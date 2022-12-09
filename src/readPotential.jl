function readPotential(potential::Potential)

    potential.internalElemEnergy = u"hartree"
    potential.internalElemCoords = u"bohr"
    potential.internalElemMass   = u"me"

    file = open(potential.file, "r")

    # removing all comments starting with #
    lines = getindex.(split.(readlines(file), "#"), 1)

    # removing all blank lines from input
    filter!(x -> !isempty(strip(x)), lines)

    lineElements = split.(lines)

    potential.dimension = length(lineElements[1]) - 1
    potential.coords    = Vector()
    potential.potential = Vector()

    [push!(potential.coords, []) for _ in 1:potential.dimension]

    for line in lineElements

        length(line)-1 != potential.dimension && (@error "Potential has inconsistent dimensionality"; exit())

        for i in 1:length(line)-1
            push!(potential.coords[i], parse(Float64, line[i]))
        end

        push!(potential.potential, parse(Float64, line[end]))

    end    

    #convert values to internal units (atomic units)
    potential.potential = ustrip.(uconvert.(potential.internalElemEnergy, potential.potential * potential.potentialUnit))
    potential.mass      = ustrip.(uconvert.(potential.internalElemMass  , potential.mass      * potential.massUnit))
    for i in eachindex(potential.coords)
        potential.coords[i]    = ustrip.(uconvert.(potential.internalElemCoords, potential.coords[i] * potential.coordsUnit))
    end

    potential.kpoints = [Tuple(zeros(potential.dimension))]

    #bandstruture setup TODO: make seperate k path routines!
    if potential.n_kpoints != -1 && potential.dimension == 1
        
        k_intervall = π / (potential.coords[1][end] - potential.coords[1][1]) / (potential.n_kpoints-1)
        [push!(potential.kpoints, k_intervall*(i-1)) for i in 1:potential.n_kpoints]

    elseif potential.n_kpoints != -1 && potential.dimension == 2

        k_intervall_1 = π / (potential.coords[1][end] - potential.coords[1][1]) / (potential.n_kpoints-1)
        k_intervall_2 = π / (potential.coords[2][end] - potential.coords[2][1]) / (potential.n_kpoints-1)

        if potential.bandStructure
            [push!(potential.kpoints, (0.0, (i-1)*k_intervall_2))                                   for i in 1:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervall_1, (potential.n_kpoints-1)*k_intervall_2)) for i in 2:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervall_1, (i-1)*k_intervall_2))                   for i in potential.n_kpoints-1:-1:1]
        else
            for i in 1:potential.n_kpoints
                kx = k_intervall_1*(i-1)
                [push!(potential.kpoints, (kx, k_intervall_2*(j-1))) for j in 1:potential.n_kpoints]
            end
        end
    elseif potential.n_kpoints != -1 && potential.dimension == 3

        k_intervall_1 = π / (potential.coords[1][end] - potential.coords[1][1]) / (potential.n_kpoints-1)
        k_intervall_2 = π / (potential.coords[2][end] - potential.coords[2][1]) / (potential.n_kpoints-1)
        k_intervall_3 = π / (potential.coords[3][end] - potential.coords[3][1]) / (potential.n_kpoints-1)

        if potential.bandStructure
            [push!(potential.kpoints, (0.0, (i-1)*k_intervall_2, 0.0)) for i in 1:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervall_1, (potential.n_kpoints-1)*k_intervall_2, 0.0)) for i in 2:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervall_1, (i-1)*k_intervall_2, 0.0)) for i in potential.n_kpoints-1:-1:1]
            [push!(potential.kpoints, ((i-1)*k_intervall_1, (i-1)*k_intervall_2, (i-1)*k_intervall_3)) for i in 2:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervall_1, (potential.n_kpoints-1)*k_intervall_2, (i-1)*k_intervall_3)) for i in potential.n_kpoints-1:-1:1]
            [push!(potential.kpoints, ((potential.n_kpoints-1)*k_intervall_1, (potential.n_kpoints-1)*k_intervall_2, (i-1)*k_intervall_3)) for i in 1:potential.n_kpoints]
        else
            for i in 1:potential.n_kpoints
                kx = k_intervall_1*(i-1)
                for j in 1:potential.n_kpoints
                    ky = k_intervall_2*(i-1)
                    [push!(potential.kpoints, (kx, ky, k_intervall_3*(k-1))) for k in 1:potential.n_kpoints]
                end
            end
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

    potential.shift     = minimum(potential.potential) #shift potential - sparse solver not capable of using not shifted potential!
    potential.potential = potential.potential .- potential.shift
    potential.intervall = potential.coords[end][2] - potential.coords[end][1]

    #check if datapoints matches input dimensions

    #check_spacing()
    #reorder ???? dont know if it should be user handled - ask thh

end