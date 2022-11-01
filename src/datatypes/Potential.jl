mutable struct Potential
    file::String

    potentialUnit::Unitful.FreeUnits
    coordsUnit   ::Unitful.FreeUnits
    massUnit     ::Unitful.FreeUnits

    internalElemEnergy::Quantity
    internalElemCoords::Quantity
    internalElemMass  ::Quantity

    dimension::Int64

    mass::Float64

    potential::Vector{Float64}

    coords::Vector{Vector{Float64}}

    Potential() = new()
end