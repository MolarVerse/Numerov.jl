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

    potential.kpoints = Vector()

    #bandstruture setup

    if potential.n_kpoints != -1 && potential.dimension == 1
        
        push!(potential.kpoints, zeros(potential.n_kpoints))
        k_intervall = π / (potential.coords[1][end] - potential.coords[1][1]) / (potential.n_kpoints-1)
        for i in 1:potential.n_kpoints
            potential.kpoints[1][i] = k_intervall*(i-1)
        end
    elseif potential.n_kpoints != -1 && potential.dimension == 2
        push!(potential.kpoints, zeros(potential.n_kpoints))
        push!(potential.kpoints, zeros(potential.n_kpoints))
        k_intervall_1 = π / (potential.coords[1][end] - potential.coords[1][1]) / (potential.n_kpoints-1)
        k_intervall_2 = π / (potential.coords[2][end] - potential.coords[2][1]) / (potential.n_kpoints-1)
        for i in 1:potential.n_kpoints
            potential.kpoints[1][i] = k_intervall_1*(i-1) #TODO: make it for arbitrary k_points
            potential.kpoints[2][i] = k_intervall_2*(i-1)
        end
    elseif potential.n_kpoints != -1 && potential.dimension == 3
        push!(potential.kpoints, zeros(potential.n_kpoints))
        push!(potential.kpoints, zeros(potential.n_kpoints))
        push!(potential.kpoints, zeros(potential.n_kpoints))
        k_intervall_1 = π / (potential.coords[1][end] - potential.coords[1][1]) / (potential.n_kpoints-1)
        k_intervall_2 = π / (potential.coords[2][end] - potential.coords[2][1]) / (potential.n_kpoints-1)
        k_intervall_3 = π / (potential.coords[3][end] - potential.coords[3][1]) / (potential.n_kpoints-1)
        for i in 1:potential.n_kpoints
            potential.kpoints[1][i] = k_intervall_1*(i-1) #TODO: make it for arbitrary k_points
            potential.kpoints[2][i] = k_intervall_2*(i-1)
            potential.kpoints[3][i] = k_intervall_3*(i-1)
        end
    else
        for i in 1:potential.dimension
            push!(potential.kpoints, zeros(1))
        end    
    end

    if isempty(potential.n_datapoints)
        if potential.dimension > 1
            @error "If the dimension is higher than one the number of datapoints per dimension has to be given explicitly in the input file e.g. datapoints = 20, 30"
            exit()
        else
            potential.n_datapoints = [length(potential.potential)]
        end
    end

    #check if datapoints matches input dimensions

    #check_spacing()
    #reorder ???? dont know if it should be user handled - ask thh

end