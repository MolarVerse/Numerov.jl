using TimerOutputs

function solve(potential::Potential, system::System, output::Output, k, to)

    @timeit to "inner solve" begin

    intervall   = potential.coords[end][2] - potential.coords[end][1]

    k_squared = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))
    Δ         = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))

    @timeit to "build k*nabla" begin
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

        k_squared = spdiagm(ones(prod(potential.n_datapoints))*(k[1]^2+k[2]^2))

    end
    end

    @timeit to "hamiltonian" Hamiltonian = 0.5 / potential.mass * (-system.laplace/intervall^2/2^(potential.dimension-1) - 2*im*Δ/intervall + k_squared) + spdiagm(potential.potential)
    
    # @timeit to "diagonalize" eigenvalues, eigenvectors = eigen(sparse(Hamiltonian))
    @timeit to "diagonalize" eigenvalues, eigenvectors = eigs(sparse(Hamiltonian), nev = output.n_eigenvalues, which = :SM)

    output.eigenvectors = Vector()
    @timeit to "assign eigvals" output.eigenvalues  = real.(eigenvalues[1:output.n_eigenvalues])
    @timeit to "assign frequencies" output.frequencies  = zeros(output.n_eigenvalues-1, output.n_eigenvalues-1)
    
    for i in 1:output.n_eigenvalues
        @timeit to "assign eigvectors" push!(output.eigenvectors, real.(eigenvectors[:,i]))
    end

    end

end