function solve(potential::Potential, system::System, output::Output, k)

    intervall   = potential.coords[end][2] - potential.coords[end][1]

    k_squared = zeros(prod(potential.n_datapoints), prod(potential.n_datapoints))
    Δ         = zeros(prod(potential.n_datapoints), prod(potential.n_datapoints))

    if potential.dimension == 1
        Δ = system.Δ*k
        k_squared = diagm(ones(system.n_datapoints)*k^2)
    elseif potential.dimension == 2
        if k[2] != 0.0
            Δ = system.Δ #* k[2]
            for i in 1:potential.n_datapoints[1]
                Δ[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] *= k[1]/k[2]
            end
        else
            Δ = system.Δ * k[2]
            Δ_tmp = copy(Δ)
            for i in 1:potential.n_datapoints[1]
                Δ[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] = Δ_tmp[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] * k[1]
            end
        end

        k_squared = diagm(ones(prod(potential.n_datapoints))*(k[1]^2+k[2]^2))

    end

    println(issymmetric(system.laplace))
    println(issymmetric(system.Δ))
    println(issymmetric(k_squared))
    println(sum(abs.(system.Δ .+ system.Δ')))

    Hamiltonian = 0.5 / potential.mass * (-system.laplace/intervall^2/2^(potential.dimension-1) - 2*im*Δ/intervall/2^(potential.dimension-1) + k_squared) + diagm(potential.potential)
    # Hamiltonian = 0.5 / potential.mass * (-system.laplace/intervall^2 - 2*im*system.Δ/intervall*k + diagm(ones(system.n_datapoints)*k^2)) + diagm(potential.potential)
    # Hamiltonian = 0.5 / potential.mass * (-system.laplace/intervall^2)/2^(potential.dimension-1) + diagm(potential.potential)

    eigenvalues, eigenvectors = eigen(Hamiltonian)

    println(eigenvalues[1:2])

    output.eigenvectors = Vector()
    output.eigenvalues  = real.(eigenvalues[1:output.n_eigenvalues])
    output.frequencies  = zeros(output.n_eigenvalues-1, output.n_eigenvalues-1)
    
    for i in 1:output.n_eigenvalues
        push!(output.eigenvectors, real.(eigenvectors[:,i]))
    end

end