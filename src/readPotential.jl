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

    potential.kpoints = []
    potential.intervall = potential.coords[end][2] - potential.coords[end][1]

    k_intervalls = [π / (potential.coords[i][end] - potential.coords[i][1] + potential.intervall) / (potential.n_kpoints-1) for i in 1:potential.dimension]

    #bandstruture setup TODO: make seperate k path routines!
    if potential.n_kpoints != -1 && potential.dimension == 1
        
        potential.kpoints = get_kpoints_1D(k_intervalls, potential.n_kpoints)

    elseif potential.n_kpoints != -1 && potential.dimension == 2

        println("test")
        potential.kpoints = get_kpoints_2D(k_intervalls, potential.n_kpoints, potential.bandStructure)
        println("test")
        
    elseif potential.n_kpoints != -1 && potential.dimension == 3

        if potential.bandStructure
            [push!(potential.kpoints, (0.0, (i-1)*k_intervalls[2], 0.0)) for i in 1:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervalls[1], (potential.n_kpoints-1)*k_intervalls[2], 0.0)) for i in 2:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervalls[1], (i-1)*k_intervalls[2], 0.0)) for i in potential.n_kpoints-1:-1:1]
            [push!(potential.kpoints, ((i-1)*k_intervalls[1], (i-1)*k_intervalls[2], (i-1)*k_intervalls[3])) for i in 2:potential.n_kpoints]
            [push!(potential.kpoints, ((i-1)*k_intervalls[1], (potential.n_kpoints-1)*k_intervalls[2], (i-1)*k_intervalls[3])) for i in potential.n_kpoints-1:-1:1]
            [push!(potential.kpoints, ((potential.n_kpoints-1)*k_intervalls[1], (potential.n_kpoints-1)*k_intervalls[2], (i-1)*k_intervalls[3])) for i in 1:potential.n_kpoints]
        else
            for i in 1:potential.n_kpoints
                kx = k_intervalls[1]*(i-1)
                for j in 1:potential.n_kpoints
                    ky = k_intervalls[2]*(i-1)
                    [push!(potential.kpoints, (kx, ky, k_intervalls[3]*(k-1))) for k in 1:potential.n_kpoints]
                end
            end
        end
    else
        push!(potential.kpoints, Tuple(zeros(potential.dimension)))
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

    #check if datapoints matches input dimensions

    #check_spacing() ???
    #reorder ???? dont know if it should be user handled - ask thh

end

function get_kpoints_1D(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    k_range = 0.0:k_intervalls[1]:k_intervalls[1]*(n_kpoints-1)

    return collect(Iterators.product(k_range))

end

function get_kpoints_2D(k_intervalls::Vector{Float64}, n_kpoints::Int64, bandStructure::Bool)

    k_points = Vector()

    if bandStructure

        Γ_X = [(0.0, k_intervalls[2]*i) for i in 0:n_kpoints-1]

        X_M = [((n_kpoints-1)*k_intervalls[1], i*k_intervalls[2]) for i in 0:n_kpoints-1]
        
        M_Γ = [(i*k_intervalls[1], i*k_intervalls[2]) for i in n_kpoints-1:-1:0]

        println(unique(vcat(Γ_X, X_M, M_Γ)))
        
        return unique(vcat(Γ_X, X_M, M_Γ))
    else

        k_range = [0:k_intervalls[i]:k_intervalls[i]*(n_kpoints - 1) for i in 1:2]

        return sort(reshape(collect(Iterators.product(k_range...)), n_kpoints^2))

    end
    
end

function get_kpoints_3D(k_intervalls::Vector{Float64}, n_kpoints::Int64, bandStructure::Bool)

    k_points = Vector()

    if bandStructure

        # [push!(potential.kpoints, ((i-1)*k_intervalls[1], (potential.n_kpoints-1)*k_intervalls[2], (i-1)*k_intervalls[3])) for i in potential.n_kpoints-1:-1:1]
        # [push!(potential.kpoints, ((potential.n_kpoints-1)*k_intervalls[1], (potential.n_kpoints-1)*k_intervalls[2], (i-1)*k_intervalls[3])) for i in 1:potential.n_kpoints]

        Γ_X = [(0.0, k_intervalls[2], 0.0) for i in 0:n_kpoints-1]

        X_M = [(i*k_intervalls[1], (n_kpoints-1)*k_intervalls[2], 0.0) for i in 0:n_kpoints-1]
        
        M_Γ = [(i*k_intervalls[1], i*k_intervalls[2], 0.0) for i in n_kpoints-1:-1:0]

        Γ_R = [(i*k_intervalls[1], i*k_intervalls[2], i*k_intervalls[3]) for i in 0:n_kpoints-1]

        #R_X = [(,(n_kpoints-1)*k_intervalls[1],) for i in n_kpoints-1:-1:0]

        #M_R

        return unique(vcat(Γ_X, X_M, M_Γ, Γ_R)) #more complicated!!!!
    else

        k_range = [0:k_intervalls[i]:k_intervalls[i]*(n_kpoints - 1) for i in 1:length(k_intervalls)]

        return sort(reshape(collect(Iterators.product(k_range...)), n_kpoints^length(k_intervalls)))

    end
    
end