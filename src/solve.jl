function solve(potential::Potential, system::System, output::Output, k, files::Files)

    ####################################################
    #                                                  #
    # setup total Hamiltonian which is then decomposed #
    #                                                  #
    ####################################################

    if system.reciprocal

        ###################################################
        #                                                 #
        # calculate k_squared matrix k² = kx² + ky² + kz² #
        #                                                 #
        ###################################################

        k_squared = spdiagm(ones(prod(potential.n_datapoints)) * (norm(k)^2))

        ###################################################
        #                                                 #
        # build ∇k matrix according to dxkx + dyky + dzkz #
        #                                                 #
        ###################################################

        ∇ = build∇_k(potential, system, k)

        @timeit files.to "build Ham" Hamiltonian = 0.5 * (-system.Δ / potential.intervall[1]^2 / 2^(potential.dimension - 1) - 2 * im * ∇ / potential.intervall[1] + k_squared) + spdiagm(potential.potential)

    else

        # non-periodic runs have k = 0 and ∇k = 0, so the Hamiltonian is real
        # symmetric - solving it in real arithmetic halves memory and work

        @timeit files.to "build Ham" Hamiltonian = 0.5 * (-system.Δ / potential.intervall[1]^2 / 2^(potential.dimension - 1)) + spdiagm(potential.potential)

    end

    #####################
    #                   #
    # solve Hamiltonian #
    #                   #
    #####################

    precond = system.solver == LOBPCG && !system.reciprocal ?
        KineticPreconditioner(potential, system) : nothing

    eigenvalues, eigenvectors = solveWrapper(system, output, files, Hamiltonian; preconditioner = precond)

    ######################################################
    #                                                    #
    # save eigenvalues and eigenvectors in output struct #
    #                                                    #
    ######################################################

    output.eigenvalues = real.(eigenvalues[1:output.n_eigenvalues])
    output.eigenvectors = Vector()
    [push!(output.eigenvectors, eigenvectors[:, i]) for i in 1:output.n_eigenvalues]

    #########################
    #                       #
    # normalizeeigenvectors #
    #                       #
    #########################

    normalize_eigenvectors(output, potential.intervall[1] / sqrt(potential.mass[1]), potential.dimension, potential) #does not work in general for different spacings in x,y,z

end

"""
Arpack shift-invert about a small negative σ. The potential is shifted so
that min(V) = 0, making the Hamiltonian positive (semi)definite: its smallest
eigenvalues are the ones closest to σ ≈ 0, so shift-invert converges in a few
iterations where the plain :SM mode needs thousands of restarts. σ sits
slightly BELOW zero so that H - σI stays safely invertible even when H itself
is exactly singular.
"""
function solve_arpack(Hamiltonian, nev::Int)
    σ = -1.0e-6 * maximum(abs, diag(Hamiltonian))
    return eigs(Hamiltonian, nev=nev, sigma=σ)
end

"""
Largest relative eigenpair residual max ‖Hx - λx‖ / (‖x‖ max(1, |λ|)).
Degenerate (near-zero) eigenvectors count as infinitely loose, so silently
collapsed solver output can never pass verification.
"""
function max_relative_residual(Hamiltonian, eigenvalues, eigenvectors, n::Int)
    r = 0.0
    for i in 1:min(n, length(eigenvalues))
        x  = view(eigenvectors, :, i)
        nx = norm(x)
        nx > sqrt(eps()) || return Inf
        r = max(r, norm(Hamiltonian * x .- eigenvalues[i] .* x) / (nx * max(1.0, abs(eigenvalues[i]))))
    end
    return r
end

function solveWrapper(system::System, output::Output, files::Files, Hamiltonian;
                      preconditioner = nothing)

    nev = output.n_eigenvalues + 5

    if system.solver == ARPACK

        @timeit files.to "Arpack" eigenvalues, eigenvectors = solve_arpack(Hamiltonian, nev)

    elseif system.solver == KRYLOV

        @timeit files.to "Krylov" eigenvalues, eigenvectors, info = eigsolve(Hamiltonian, nev, :SR; ishermitian=true, maxiter=10000)
        # stack the eigenvectors as matrix columns; the previous adjoint-based
        # reshape conjugated complex eigenvectors
        eigenvectors = stack(eigenvectors)

    elseif system.solver == LOBPCG

        eltype(Hamiltonian) <: Real ||
            throw(ArgumentError("the lobpcg solver supports non-periodic (real symmetric) problems - use arpack or krylov for periodic k-point runs"))

        # LOBPCG occasionally breaks down (its internal factorizations fail on
        # ill-conditioned iteration blocks), so retry with a fresh random
        # block and fall back to Arpack shift-invert if it keeps failing -
        # accuracy and robustness are never worse than the arpack path
        result = nothing
        @timeit files.to "LOBPCG" for attempt in 1:2
            try
                X0 = randn(size(Hamiltonian, 1), nev)
                result = preconditioner === nothing ?
                    lobpcg(Hamiltonian, false, X0; tol = 1.0e-7, maxiter = 2000) :
                    lobpcg(Hamiltonian, false, X0; P = preconditioner, tol = 1.0e-7, maxiter = 2000)
                break
            catch err
                @warn "lobpcg attempt $attempt failed, $(attempt == 1 ? "retrying" : "falling back to arpack")" err
                result = nothing
            end
        end

        if result === nothing
            @timeit files.to "Arpack" eigenvalues, eigenvectors = solve_arpack(Hamiltonian, nev)
        else
            eigenvalues, eigenvectors = result.λ, result.X
        end

    elseif system.solver == GPU

        throw(ArgumentError("cuda solver is not implemented"))

    else

        @timeit files.to "LU" eigenvalues, eigenvectors = eigen(Matrix(Hamiltonian))

    end

    # iterative solvers do not guarantee an ordering - sort ascending by
    # real part so downstream truncation always keeps the lowest states
    perm = sortperm(eigenvalues; by = real)
    eigenvalues, eigenvectors = eigenvalues[perm], eigenvectors[:, perm]

    # verify, don't trust: check the residuals of the eigenpairs that will be
    # kept, and escalate or warn if an iterative solver returned loose pairs
    if system.solver in (ARPACK, KRYLOV, LOBPCG)
        residual = max_relative_residual(Hamiltonian, eigenvalues, eigenvectors, output.n_eigenvalues)
        if system.solver == LOBPCG && residual > 1.0e-6
            @warn "lobpcg eigenpairs exceed the residual tolerance, re-solving with arpack" residual
            eigenvalues, eigenvectors = solve_arpack(Hamiltonian, nev)
            perm = sortperm(eigenvalues; by = real)
            eigenvalues, eigenvectors = eigenvalues[perm], eigenvectors[:, perm]
        elseif residual > 1.0e-6
            @warn "eigenpair residuals are larger than expected" solver = system.solver residual
        end
    end

    return eigenvalues, eigenvectors
end

