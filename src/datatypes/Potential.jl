mutable struct Potential
    file::String

    shift::Bool

    potentialUnit::Unitful.FreeUnits
    coordsUnit   ::Unitful.FreeUnits
    massUnit     ::Unitful.FreeUnits

    internalElemEnergy::Unitful.FreeUnits
    internalElemCoords::Unitful.FreeUnits
    internalElemMass  ::Unitful.FreeUnits

    n_kpoints::Int64
    kpoints  ::Vector{Vector{Float64}}

    dimension::Int64

    mass::Float64

    potential::Vector{Float64}

    coords::Vector{Vector{Float64}}

    n_datapoints::Vector{Int64}

    Potential() = new()
end