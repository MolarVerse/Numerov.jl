function solve(potential::Potential, system::System, output::Output, k, files::Files)

    intervall    = potential.intervall

    k_squared = spdiagm(ones(prod(potential.n_datapoints))*(norm(k)^2))

    ∇ = build∇_k(potential, system, k)

    @timeit files.to "build Ham" Hamiltonian = 0.5 / potential.mass[1] * (-system.Δ/intervall^2/2^(potential.dimension-1) - 2*im*∇/intervall + k_squared) + spdiagm(potential.potential)
    
    output.eigenvectors = Vector()

    eigenvalues, eigenvectors = solveWrapper(system, output, files, Hamiltonian)

    output.eigenvalues  = real.(eigenvalues[1:output.n_eigenvalues])
    output.frequencies  = zeros(output.n_eigenvalues-1, output.n_eigenvalues-1)
        
    [push!(output.eigenvectors, eigenvectors[:,i]) for i in 1:output.n_eigenvalues]

    normalize_eigenvectors(output, intervall, potential.dimension, potential)

end

function solveWrapper(system::System, output::Output, files::Files, Hamiltonian)

    if system.solver == ARPACK
        if isempty(output.eigenvectors)
            @timeit files.to "Arpack" eigenvalues, eigenvectors = eigs(sparse(Hamiltonian), nev = output.n_eigenvalues+5, which = :SM, maxiter=typemax(Int))
        else
            @timeit files.to "Arpack" eigenvalues, eigenvectors = eigs(sparse(Hamiltonian), nev = output.n_eigenvalues+5, which = :SM, maxiter=typemax(Int), v0=output.eigenvectors[1])
        end
    elseif system.solver == KRYLOV
        #check convergence of all eigenvalues with info
        @timeit files.to "Krylov" eigenvalues, eigenvectors, info = eigsolve(sparse(Hamiltonian), output.n_eigenvalues+5, :SR; ishermitian=true, maxiter=10000)
        @show(info)
        eigenvectors = mapreduce(permutedims, vcat, eigenvectors)'
    elseif system.solver == GPU
        @error "CUDA sparse solver not working yet!"
        exit()

        # @timeit files.to "diagonalize CUDA" begin
        #     # CUDA.device!(1)
        #     T = ComplexF32
        #     Hamiltonian = SparseMatrixCSC{ComplexF32, Int32}(Hamiltonian)
        #     Hamiltonian = CUSPARSE.CuSparseMatrixCSR(Hamiltonian)
        #     CUSOLVER.csreigvsi(Hamiltonian, rand(T), CUDA.rand(T, prod(potential.n_datapoints)), Float32(1e-6), Cint(1000), 'O')
        # end

    end

    return eigenvalues, eigenvectors
end

function normalize_eigenvectors(output::Output, intervall, dimension, potential)
    for i in 1:output.n_eigenvalues
        density = output.eigenvectors[i].^2
    
        norm = sum(density) * ustrip(uconvert(potential.coordsUnit, intervall*potential.internalElemCoords))^dimension

        output.eigenvectors[i] = output.eigenvectors[i] ./ sqrt(norm)
    end
end