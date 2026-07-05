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

    Potential() = new("",                                 # file
                      false, false, false,                # reciprocal, bandStructure, read_kpoints
                      [false],                            # periodic
                      u"hartree", u"angstrom", u"u",      # potentialUnit, coordsUnit, massUnit
                      u"hartree", u"bohr", u"m_e",         # internalElemEnergy, internalElemCoords, internalElemMass
                      -1, 0,                              # n_kpoints, dimension
                      0.0,                                # shift
                      Vector{Int64}(),                    # n_datapoints
                      [1.0],                              # mass
                      Vector{Float64}(),                  # intervall
                      Vector{Float64}(),                  # potential
                      Vector{Tuple}(),                    # kpoints
                      Vector{Vector{Float64}}())          # coords
end
