function readPotential(potential::Potential)

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
            push!(potential.coords[i], parse(Float64, line[i])*potential.coordsUnit)
        end

        push!(potential.potential, parse(Float64, line[end])*potential.potentialUnit)

    end    
end