mutable struct Output
    n_eigenvalues::Int64
    
    eigenvalues::Vector{Float64}

    eigenvectors::Vector{Vector{ComplexF64}}

    frequencies::Matrix{Float64}

    Output() = new()
end