function solve(potential::Potential, system::System, output::Output, k, files::Files)

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

    ####################################################
    #                                                  #
    # setup total Hamiltonian which is then decomposed #
    #                                                  #
    ####################################################

    @timeit files.to "build Ham" Hamiltonian = 0.5 * (-system.Δ / potential.intervall[1]^2 / 2^(potential.dimension - 1) - 2 * im * ∇ / potential.intervall[1] + k_squared) + spdiagm(potential.potential)

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

        @timeit files.to "Arpack" eigenvalues, eigenvectors = eigs(sparse(Hamiltonian), nev=output.n_eigenvalues + 5, which=:SM, maxiter=typemax(Int))

    elseif system.solver == KRYLOV

        @timeit files.to "Krylov" eigenvalues, eigenvectors, info = eigsolve(sparse(Hamiltonian), output.n_eigenvalues + 5, :SR; ishermitian=true, maxiter=10000)
        eigenvectors = mapreduce(permutedims, vcat, eigenvectors)'

    elseif system.solver == GPU

        throw(ArgumentError("cuda solver is not implemented"))

    else

        @timeit files.to "LU" eigenvalues, eigenvectors = eigen(Matrix(Hamiltonian))

    end

    return eigenvalues, eigenvectors
end

