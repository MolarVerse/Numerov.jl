using TimerOutputs
using KrylovKit
using CUDA
using CUDA:CUSPARSE

function solve(potential::Potential, system::System, output::Output, k, to, files::Files)

    @timeit to "inner solve" begin

    intervall   = potential.coords[end][2] - potential.coords[end][1]

    if system.reciprocal
        k_squared, Δ = build_hamiltonian_components(potential, system, k)
    else
        k_squared = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))
        Δ         = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))
    end

    @timeit to "hamiltonian" Hamiltonian = 0.5 / potential.mass * (-system.laplace/intervall^2/2^(potential.dimension-1) - 2*im*Δ/intervall + k_squared) + spdiagm(potential.potential)
    
    if sum(k) == 0.0
        println(files.logFile, "Non zeros = ", length(Hamiltonian.nzval))
        println(files.logFile, "zeros     = ", length(Hamiltonian) - length(Hamiltonian.nzval))
        println(files.logFile, "Sparsity  = ", length(Hamiltonian.nzval) / length(Hamiltonian))
    end

    eigenvalues, eigenvectors = solveWrapper(system, output, files, Hamiltonian, to)

    

    output.eigenvectors = Vector()
    @timeit to "assign eigvals" output.eigenvalues  = real.(eigenvalues[1:output.n_eigenvalues])
    @timeit to "assign frequencies" output.frequencies  = zeros(output.n_eigenvalues-1, output.n_eigenvalues-1)
    
    for i in 1:output.n_eigenvalues
        @timeit to "assign eigvectors" push!(output.eigenvectors, real.(eigenvectors[:,i]))
    end

    normalize_eigenvectors(output, intervall, potential.dimension, potential)

    end

end

function solveWrapper(system::System, output::Output, files::Files, Hamiltonian, to)

    if system.solver == ARPACK
        @timeit to "Arpack" eigenvalues, eigenvectors = eigs(sparse(Hamiltonian), nev = output.n_eigenvalues+5, which = :SM, maxiter=typemax(Int))
    elseif system.solver == KRYLOV
        @timeit to "Krylov" eigenvalues, eigenvectors, info = eigsolve(sparse(Hamiltonian), output.n_eigenvalues+5, :SR; ishermitian=true, maxiter=10000)
        @show(info)
        eigenvectors = mapreduce(permutedims, vcat, eigenvectors)'
    elseif system.solver == GPU
        @error "CUDA sparse solver not working yet!"
        exit()

        # @timeit to "diagonalize CUDA" begin
        #     # CUDA.device!(1)
        #     T = ComplexF32
        #     Hamiltonian = SparseMatrixCSC{ComplexF32, Int32}(Hamiltonian)
        #     Hamiltonian = CUSPARSE.CuSparseMatrixCSR(Hamiltonian)
        #     CUSOLVER.csreigvsi(Hamiltonian, rand(T), CUDA.rand(T, prod(potential.n_datapoints)), Float32(1e-6), Cint(1000), 'O')
        # end

    end

    return eigenvalues, eigenvectors
end

using NumericalIntegration

function normalize_eigenvectors(output::Output, intervall, dimension, potential)
    for i in 1:output.n_eigenvalues
        density = output.eigenvectors[i].^2
    
        norm = sum(density) * ustrip(uconvert(potential.coordsUnit, intervall*potential.internalElemCoords))^dimension

        output.eigenvectors[i] = output.eigenvectors[i] ./ sqrt(norm)
    end
end

function build_hamiltonian_components(potential::Potential, system::System, k)
    k_squared = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))
    Δ         = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))

    if potential.dimension == 1
        Δ = system.Δ*k
        k_squared = spdiagm(ones(system.n_datapoints)*k^2)
    elseif potential.dimension == 2
        if k[2] != 0.0
            Δ = system.Δ * k[2]
            for i in 1:potential.n_datapoints[1]
                Δ[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] *= k[1]/k[2]
            end
        else
            for i in 1:potential.n_datapoints[1]
                Δ[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] = system.Δ[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] * k[1]
            end
        end

        k_squared = spdiagm(ones(prod(potential.n_datapoints))*(norm(k)^2))

    elseif potential.dimension == 3

        stencil    = zeros(system.stencil, system.stencil, system.stencil)

        stencil[:                 ,system.stencil÷2+1,system.stencil÷2+1] = ones(system.stencil)*k[1]
        stencil[system.stencil÷2+1,:                 ,system.stencil÷2+1] = ones(system.stencil)*k[2]
        stencil[system.stencil÷2+1,system.stencil÷2+1,:                 ] = ones(system.stencil)*k[3]
        stencil[system.stencil÷2+1,system.stencil÷2+1,system.stencil÷2+1] = 0.0
 
        Δ = system.Δ .* build_3d_stencil(system, potential.n_datapoints, stencil)

        k_squared = spdiagm(ones(prod(potential.n_datapoints))*(norm(k)^2))

    end

    return k_squared, Δ
end
