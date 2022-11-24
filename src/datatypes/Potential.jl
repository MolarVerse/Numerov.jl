mutable struct Potential
    file::String

    shift        ::Bool
    reciprocal   ::Bool
    bandStructure::Bool

    potentialUnit::Unitful.FreeUnits
    coordsUnit   ::Unitful.FreeUnits
    massUnit     ::Unitful.FreeUnits

    internalElemEnergy::Unitful.FreeUnits
    internalElemCoords::Unitful.FreeUnits
    internalElemMass  ::Unitful.FreeUnits

    n_kpoints::Int64
    kpoints  ::Vector{}

    dimension::Int64

    mass::Float64

    potential::Vector{Float64}

    coords::Vector{Vector{Float64}}

    n_datapoints::Vector{Int64}

    Potential() = new()
end