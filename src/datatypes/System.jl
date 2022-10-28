mutable struct System1D <: System
    stencil     ::Int64
    n_datapoints::Int64

    laplace::SparseMatrixCSC{Float64, Int64}

    Δ::Matrix{Quantity}

    System1D() = new()
end