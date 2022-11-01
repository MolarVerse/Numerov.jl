function readPotential(potential::Potential)

    potential.internalElemEnergy = u"hartree"
    potential.internalElemCoords = u"bohr"
    potential.internalElemMass   = u"me"

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

    potential.potential = ustrip.(uconvert.(potential.internalElemEnergy, potential.potential * potential.potentialUnit))
    potential.mass      = ustrip.(uconvert.(potential.internalElemMass  , potential.mass      * potential.massUnit)) #still no idea why this horrible hack with .)

    for i in eachindex(potential.coords)
        potential.coords[i]    = ustrip.(uconvert.(potential.internalElemCoords, potential.coords[i] * potential.coordsUnit))
    end

    #shift potential

    potential.shift && (potential.potential = potential.potential .- minimum(potential.potential))
end