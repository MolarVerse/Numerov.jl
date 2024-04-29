function readPotential(potential::Potential, files::Files)

    #######################################
    #                                     #
    # setup internal units (atomic units) #
    #                                     #
    #######################################

    set_internalUnits(potential)

    ######################################################
    #                                                    #
    # prepare potential file for processing to read data #
    #                                                    #
    ######################################################

    file = open(files.PotentialFileName, "r")

    lines = getindex.(split.(readlines(file), "#"), 1) # removing all comments starting with #
    
    filter!(x -> !isempty(strip(x)), lines) # removing all blank lines from input

    lineElements = split.(lines)

    ############################################################
    #                                                          #
    # determine dimension of problem and initalize all vectors #
    #                                                          #
    ############################################################

    potential.dimension = length(lineElements[1]) - 1
    potential.coords    = Vector()
    potential.potential = Vector()

    [push!(potential.coords, []) for _ in 1:potential.dimension]

    ##############################################################
    #                                                            #
    # read lines from potential file to get potential and coords #
    #                                                            #
    ##############################################################

    for line in lineElements

        length(line)-1 != potential.dimension && (@error "Potential has inconsistent dimensionality"; exit())

        for i in 1:length(line)-1
            push!(potential.coords[i], parse(Float64, line[i]))
        end

        push!(potential.potential, parse(Float64, line[end]))

    end    

    ###############################################################################################################
    #                                                                                                             #
    # repeat periodic for higher dimensions - if only one periodic keyword given -> same periodicity for all DOFs #
    #                                                                                                             #
    # if more than one periodic keyword is given - periodicity for all dimensions has to be given                 #
    #                                                                                                             #
    ###############################################################################################################

    length(potential.periodic) == 1 && (potential.periodic = repeat(potential.periodic, potential.dimension))

    length(potential.periodic) != potential.dimension && (@error "periodic is not correctly defined regarding the potential dimension!"; exit())

    ##########################################################################################
    #                                                                                        #
    # repeat masses for higher dimensions - if only one mass given -> same mass for all DOFs #
    #                                                                                        #
    # if more than one mass is given - masses for all dimensions have to be given            #
    #                                                                                        #
    ##########################################################################################

    length(potential.mass) == 1 && (potential.mass = repeat(potential.mass, potential.dimension))

    length(potential.mass) != potential.dimension && (@error "reduced-mass is not correctly defined regarding the potential dimension!"; exit())

    #######################################################################
    #                                                                     #
    # convert potential, coords and mass to internal units (atomic units) #
    #                                                                     #
    #######################################################################

    convert_to_internalUnits(potential)

    #################################################
    #                                               #
    # determine all unique coordinates for all DOFs # TODO: check spacing with tolerance keyword!
    #                                               #
    #################################################

    uniquecoords = []
    [push!(uniquecoords, unique(potential.coords[i])) for i in 1:potential.dimension]

    ###################################################
    #                                                 #
    # calculate mass weighted intervalls for all DOFs #
    #                                                 #
    ###################################################

    potential.intervall = zeros(potential.dimension)
    [potential.intervall[i] = (uniquecoords[i][2] - uniquecoords[i][1]) * sqrt(potential.mass[i]) for i in 1:potential.dimension]

    #############TODO:###############

    if files.k_pointsFileName != "" && !potential.read_kpoints
        @warn "You have defined a k-points file but not set read kpoints to true"
    end #TODO: write reading logic and integrate it

    ########################################
    #                                      #
    # calculate mass weighted k intervalls #
    #                                      #
    ########################################

    k_intervalls = [π / ((potential.coords[i][end] - potential.coords[i][1]) * sqrt(potential.mass[i]) + potential.intervall[i]) / (potential.n_kpoints-1) for i in 1:potential.dimension]

    potential.n_kpoints != -1 ? potential.reciprocal = true : potential.reciprocal = false

    #############################################################################################################################
    #                                                                                                                           #
    # determine all k-points for which Schrödinger equation has to be solved                                                    #
    #                                                                                                                           #
    # depending on input no k-points are included or all possible combinations or the minimal k-path through the brioullin zone #
    #                                                                                                                           #
    #############################################################################################################################

    if potential.n_kpoints != -1 
        if potential.bandStructure
            if potential.dimension == 1
        
                potential.kpoints = get_kpoints_1D(k_intervalls, potential.n_kpoints)

            elseif potential.dimension == 2

                potential.kpoints = get_kpoints_2D(k_intervalls, potential.n_kpoints)
        
            elseif potential.dimension == 3

                potential.kpoints = get_kpoints_3D(k_intervalls, potential.n_kpoints)

            end
        else

            k_range = [0:k_intervalls[i]:k_intervalls[i]*(potential.n_kpoints - 1) for i in eachindex(k_intervalls)]

            potential.kpoints = sort(reshape(collect(Iterators.product(k_range...)), potential.n_kpoints^length(k_intervalls)))
        end
    else
        potential.kpoints = [Tuple(zeros(potential.dimension))]
    end

    #######################################################
    #                                                     #
    # filter for selecting kpoints only for periodic DOFs #
    #                                                     #
    #######################################################

    for i in 1:potential.dimension
        if !potential.periodic[i]
            potential.kpoints = filter(x -> x[i] == 0.0, potential.kpoints)
        end
    end

    ################################################################################################################
    #                                                                                                              #
    # setup datapoints for calculation - for higher dimensions n_datapoints has to be explicitly set in input file #
    #                                                                                                              #
    ################################################################################################################

    if isempty(potential.n_datapoints)
        if potential.dimension > 1
            @error "If the dimension is higher than one the number of datapoints per dimension has to be given explicitly in the input file e.g. datapoints = 20, 30"
            exit()
        else
            potential.n_datapoints = [length(potential.potential)]
        end
    end

    ####################################################################
    #                                                                  #
    # pot_min for shifting potential later - necessary for eigs solver #
    #                                                                  #
    ####################################################################

    potential.shift = minimum(potential.potential)

end

function set_internalUnits(potential::Potential)
    potential.internalElemEnergy = u"hartree"
    potential.internalElemCoords = u"bohr"
    potential.internalElemMass   = u"me"
end

function convert_to_internalUnits(potential::Potential)
    potential.potential = ustrip.(uconvert.(potential.internalElemEnergy, potential.potential * potential.potentialUnit))
    potential.mass      = ustrip.(uconvert.(potential.internalElemMass  , potential.mass      * potential.massUnit))
    for i in eachindex(potential.coords)
        potential.coords[i]    = ustrip.(uconvert.(potential.internalElemCoords, potential.coords[i] * potential.coordsUnit))
    end
end