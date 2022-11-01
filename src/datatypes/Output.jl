mutable struct Output
    n_eigenvalues::Int64
    
    eigenvalues::Vector{Float64}

    eigenvectors::Vector{Vector{Float64}}

    frequencies::Matrix{Float64}

    Output() = new()
end