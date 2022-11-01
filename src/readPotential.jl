function readPotential(potential::Potential)

    potential.internalElemEnergy = 1.0u"hartree"
    potential.internalElemCoords = 1.0u"bohr"
    potential.internalElemMass   = 1.0u"me"

    file = open(potential.file, "r")

    lines = readlines(file)

    # removing all comments starting with #
    lines = getindex.(split.(lines, "#"), 1)

    # removing all blank lines from input
    filter!(x -> !isempty(strip(x)), lines)

    lineElements = split.(lines)

    potential.dimension = length(lineElements[1]) - 1

    potential.coords    = Vector()
    potential.potential = Vector()

    for _ in 1:potential.dimension
        push!(potential.coords, [])
    end

    for line in lineElements

        length(line)-1 != potential.dimension && (@error "Potential has inconsistent dimensionality"; exit())

        for i in 1:length(line)-1
            push!(potential.coords[i], parse(Float64, line[i]))
        end

        push!(potential.potential, parse(Float64, line[end]))

    end    

    potential.potential = ustrip.(uconvert.(unit(potential.internalElemEnergy), potential.potential * potential.potentialUnit))
    potential.mass      = ustrip.(uconvert.(unit(potential.internalElemMass  ), potential.mass      * potential.massUnit)) #still no idea why this horrible hack with unit(.)

    for i in eachindex(potential.coords)
        potential.coords[i]    = ustrip.(uconvert.(unit(potential.internalElemCoords), potential.coords[i]    * potential.coordsUnit))
    end
end