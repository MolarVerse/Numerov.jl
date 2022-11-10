function solve(potential::Potential, system::System, output::Output, k::Float64)

    intervall   = potential.coords[1][2] - potential.coords[1][1]

    Hamiltonian = 0.5 / potential.mass * (-system.laplace/intervall^2 - 2*im*system.Δ/intervall*k + diagm(ones(system.n_datapoints)*k^2)) + diagm(potential.potential)

    eigenvalues, eigenvectors = eigen(Hamiltonian)

    output.eigenvectors = Vector()
    output.eigenvalues  = eigenvalues[1:output.n_eigenvalues]
    output.frequencies  = zeros(output.n_eigenvalues-1, output.n_eigenvalues-1)
    
    for i in 1:output.n_eigenvalues
        push!(output.eigenvectors, real.(eigenvectors[:,i]))
    end

end