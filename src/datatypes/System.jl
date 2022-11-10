mutable struct System1D <: System
    stencil     ::Int64
    n_datapoints::Int64

    periodic     ::Bool
    bandStructure::Bool

    laplace::Matrix{Float64}

    Δ::Matrix{Float64}

    System1D() = new()
end

mutable struct System2D <: System
    stencil     ::Int64
    n_datapoints::Vector{Int64}

    periodic     ::Bool
    bandStructure::Bool

    laplace::Matrix{Float64}

    Δ::Matrix{Float64}

    System2D() = new()
end