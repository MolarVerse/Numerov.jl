mutable struct Potential
    file::String

    periodic     ::Vector{Bool}
    reciprocal   ::Bool
    bandStructure::Bool

    potentialUnit::Unitful.FreeUnits
    coordsUnit   ::Unitful.FreeUnits
    massUnit     ::Unitful.FreeUnits

    internalElemEnergy::Unitful.FreeUnits
    internalElemCoords::Unitful.FreeUnits
    internalElemMass  ::Unitful.FreeUnits

    n_kpoints::Int64
    kpoints  ::Vector{Tuple}

    dimension::Int64

    mass::Vector{Float64}

    shift::Float64

    intervall::Vector{Float64}

    potential::Vector{Float64}

    coords::Vector{Vector{Float64}}

    n_datapoints::Vector{Int64}

    Potential() = new()
end