mutable struct System
    stencil     ::Int64
    stencil∇    ::Int64
    stencilΔ    ::Int64
    n_datapoints::Vector{Int64}

    periodic     ::Vector{Bool}
    reciprocal   ::Bool
    bandStructure::Bool

    ∇::SparseMatrixCSC{Float64, Int64}
    Δ::SparseMatrixCSC{Float64, Int64}

    solver::SolverEnum

    System() = new()
end