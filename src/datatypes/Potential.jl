mutable struct Potential
    file::String

    reciprocal   ::Bool
    bandStructure::Bool
    read_kpoints ::Bool

    periodic     ::Vector{Bool}

    potentialUnit::Unitful.FreeUnits
    coordsUnit   ::Unitful.FreeUnits
    massUnit     ::Unitful.FreeUnits

    internalElemEnergy::Unitful.FreeUnits
    internalElemCoords::Unitful.FreeUnits
    internalElemMass  ::Unitful.FreeUnits

    n_kpoints::Int64
    dimension::Int64

    shift::Float64

    n_datapoints::Vector{Int64}
    mass        ::Vector{Float64}
    intervall   ::Vector{Float64}
    potential   ::Vector{Float64}
    kpoints     ::Vector{Tuple}

    coords::Vector{Vector{Float64}}

    Potential() = new()
end