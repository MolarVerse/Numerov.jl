mutable struct System1D <: System
    stencil     ::Int64
    n_datapoints::Int64

    periodic     ::Vector{Bool}
    reciprocal   ::Bool
    bandStructure::Bool

    laplace::SparseMatrixCSC{Float64, Int64}

    Δ::SparseMatrixCSC{Float64, Int64}

    System1D() = new()
end

mutable struct System2D <: System
    stencil     ::Int64
    n_datapoints::Vector{Int64}

    periodic     ::Vector{Bool}
    reciprocal   ::Bool
    bandStructure::Bool

    laplace::SparseMatrixCSC{Float64, Int64}

    Δ::SparseMatrixCSC{Float64, Int64}

    System2D() = new()
end

mutable struct System3D <: System
    stencil     ::Int64
    n_datapoints::Vector{Int64}

    periodic     ::Vector{Bool}
    reciprocal   ::Bool
    bandStructure::Bool

    laplace::SparseMatrixCSC{Float64, Int64}

    Δ::SparseMatrixCSC{Float64, Int64}

    System3D() = new()
end