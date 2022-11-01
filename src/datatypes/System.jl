mutable struct System1D <: System
    stencil     ::Int64
    n_datapoints::Int64

    periodic ::Bool

    laplace::SparseMatrixCSC{Float64, Int64}

    Δ::Matrix{Float64}

    System1D() = new()
end