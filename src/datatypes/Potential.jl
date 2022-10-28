mutable struct Potential
    file::String

    potentialUnit::Unitful.FreeUnits
    coordsUnit   ::Unitful.FreeUnits

    dimension::Int64
    periodic ::Bool

    potential::Vector{Quantity}

    coords::Vector{Vector{Quantity}}

    Potential() = new()
end