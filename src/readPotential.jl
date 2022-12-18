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

    uniquecoords = unique(potential.coords[i] for i in 1:potential.dimension)

    potential.intervall = zeros(potential.dimension)

    println(uniquecoords)
    
    [potential.intervall[i] = uniquecoords[i][2] - uniquecoords[i][1] for i in 1:potential.dimension]

    k_intervalls = [π / (potential.coords[i][end] - potential.coords[i][1] + potential.intervall[i]) / (potential.n_kpoints-1) for i in 1:potential.dimension]

    if length(potential.periodic) == 1
        potential.periodic = repeat(potential.periodic, potential.dimension)
    end

    #mass check
    if length(potential.mass) == 1
        potential.mass = repeat(potential.mass, potential.dimension)
    end

    if length(potential.mass) != potential.dimension
        @error "reduced-mass is not correctly defined regarding the potential dimension!"
        exit()
    end

    potential.n_kpoints != -1 ? potential.reciprocal = true : potential.reciprocal = false

    #bandstruture setup TODO: make seperate k path routines!
    if potential.n_kpoints != -1 
        if potential.bandStructure
            if potential.dimension == 1
        
                potential.kpoints = get_kpoints_1D(k_intervalls, potential.n_kpoints)

            elseif potential.dimension == 2

                potential.kpoints = get_kpoints_2D(k_intervalls, potential.n_kpoints, potential.bandStructure)
        
            elseif potential.dimension == 3

                potential.kpoints = get_kpoints_3D(k_intervalls, potential.n_kpoints, potential.bandStructure)

            end
        else

            k_range = [0:k_intervalls[i]:k_intervalls[i]*(potential.n_kpoints - 1) for i in eachindex(k_intervalls)]

            potential.kpoints = sort(reshape(collect(Iterators.product(k_range...)), potential.n_kpoints^length(k_intervalls)))
        end
    else
        potential.kpoints = [Tuple(zeros(potential.dimension))]
    end

    #filter only periodic directions of potential!
    for i in 1:potential.dimension
        if !potential.periodic[i]
            potential.kpoints = filter(x -> x[i] == 0.0, potential.kpoints)
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

    #check if datapoints matches input dimensions

    #check_spacing() ???
    #reorder ???? dont know if it should be user handled - ask thh

end

function get_kpoints_1D(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    k_range = 0.0:k_intervalls[1]:k_intervalls[1]*(n_kpoints-1)

    return collect(Iterators.product(k_range))

end

function get_kpoints_2D(k_intervalls::Vector{Float64}, n_kpoints::Int64, bandStructure::Bool)

    Γ_X = [(k_intervalls[1]*i, 0.0) for i in 0:n_kpoints-1]

    X_M = [((n_kpoints-1)*k_intervalls[1], i*k_intervalls[2]) for i in 0:n_kpoints-1]
    
    M_Γ = [(i*k_intervalls[1], i*k_intervalls[2]) for i in n_kpoints-1:-1:0]

    return rle(vcat(Γ_X, X_M, M_Γ))[1]
end

function get_kpoints_3D(k_intervalls::Vector{Float64}, n_kpoints::Int64, bandStructure::Bool)

    Γ_X = [(k_intervalls[1]*i, 0.0, 0.0) for i in 0:n_kpoints-1]

    X_M = [((n_kpoints-1)*k_intervalls[1], i*k_intervalls[2], 0.0) for i in 0:n_kpoints-1]
    
    M_Γ = [(i*k_intervalls[1], i*k_intervalls[2], 0.0) for i in n_kpoints-1:-1:0]

    Γ_R = [(i*k_intervalls[1], i*k_intervalls[2], i*k_intervalls[3]) for i in 0:n_kpoints-1]

    R_X = [((n_kpoints-1)*k_intervalls[1], i*k_intervalls[2], i*k_intervalls[3]) for i in n_kpoints-1:-1:0]
    
    M_R = [(0.0, 0.0, k_intervalls[3]*i) for i in 0:n_kpoints-1]
    
    return rle(vcat(Γ_X, X_M, M_Γ, Γ_R, R_X, M_R))[1]
end