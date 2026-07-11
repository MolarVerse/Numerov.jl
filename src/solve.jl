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

    eigenvalues, eigenvectors = solveWrapper(system, output, files, Hamiltonian)

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

function solveWrapper(system::System, output::Output, files::Files, Hamiltonian)

    if system.solver == ARPACK

        # The potential is shifted so that min(V) = 0, making the Hamiltonian
        # positive (semi)definite: its smallest eigenvalues are the ones
        # closest to σ ≈ 0, so shift-invert mode converges in a few
        # iterations where the plain :SM mode needs thousands of restarts.
        # σ is placed slightly BELOW zero so that H - σI stays safely
        # invertible even when H itself is exactly singular.
        σ = -1.0e-6 * maximum(abs, diag(Hamiltonian))
        @timeit files.to "Arpack" eigenvalues, eigenvectors = eigs(Hamiltonian, nev=output.n_eigenvalues + 5, sigma=σ)

    elseif system.solver == KRYLOV

        @timeit files.to "Krylov" eigenvalues, eigenvectors, info = eigsolve(Hamiltonian, output.n_eigenvalues + 5, :SR; ishermitian=true, maxiter=10000)
        # stack the eigenvectors as matrix columns; the previous adjoint-based
        # reshape conjugated complex eigenvectors
        eigenvectors = stack(eigenvectors)

    elseif system.solver == GPU

        throw(ArgumentError("cuda solver is not implemented"))

    else

        @timeit files.to "LU" eigenvalues, eigenvectors = eigen(Matrix(Hamiltonian))

    end

    # iterative solvers do not guarantee an ordering - sort ascending by
    # real part so downstream truncation always keeps the lowest states
    perm = sortperm(eigenvalues; by = real)

    return eigenvalues[perm], eigenvectors[:, perm]
end

