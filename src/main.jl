function numerov(inputFileName::String)

    #########################################################################
    #                                                                       #
    # initialize structs and timeroutput + reset inputdictionary to default #
    #                                                                       #
    #########################################################################

    potential = Potential()
    system    = System1D() #default setup but gets overridden later!
    output    = Output()
    files     = Files()

    files.to = TimerOutput()

    [inputDictionary[key] = "" for key in keys(inputDictionary)] #to reset dict to default values if calculation started in same repl session

    @timeit files.to "main" begin

        ###################
        #                 #
        # read input file #
        #                 #
        ###################

        readInputFile(inputFileName)
        
        #########################
        #                       #
        # parse and check input #
        #                       #
        #########################

        checkInput(potential)
        checkInput(system)
        checkInput(files)
        checkInput(output)

        #######################
        #                     #
        # read potential file #
        #                     #
        #######################

        readPotential(potential)

        ################
        #              #
        # setup system #
        #              #
        ################

        system = setupSystem(potential, system)

        ################################################################################
        #                                                                              #
        # build ∇ and Δ matrix - caution ∇ matrix later modified according to k-points #
        #                                                                              #
        ################################################################################

        buildΔ(system)
        build∇(system)

        ##########################################
        #                                        #
        # print sparsity information to log file #
        #                                        #
        ##########################################

        println(files.logFile, "Non zeros = ", length(system.Δ.nzval))
        println(files.logFile, "zeros     = ", length(system.Δ) - length(system.Δ.nzval))
        println(files.logFile, "Sparsity  = ", length(system.Δ.nzval) / length(system.Δ))

        isfile("eigenvalues.dat") && rm("eigenvalues.dat") #rm eigenvalue file if it exists TODO: think of a way to restart calculation for different k

        ########################################################################################
        #                                                                                      #
        # loop over all k-points (for non reciprocal system only one single point calculation) #
        #                                                                                      #
        ########################################################################################

        @timeit files.to "loop" begin
            for (i, k) in enumerate(potential.kpoints)

                #########################################
                #                                       #
                # shift potential to pot_min equals 0.0 #
                #                                       #
                #########################################

                potential.potential = potential.potential .- potential.shift

                ##############################
                #                            #
                # solve Schrödinger equation #
                #                            #
                ##############################

                @timeit files.to "solve" solve(potential, system, output, k, files)

                ####################################
                #                                  #
                # calculate k from mass weighted k #
                #                                  #
                ####################################

                k = k .* sqrt.(potential.mass)

                #############################################################
                #                                                           #
                # setup all file names depending on the momentanous k-point #
                #                                                           #
                #############################################################

                k_string = join(ustrip.(uconvert.(potential.coordsUnit^(-1), k ./ potential.internalElemCoords)), "_")

                if system.reciprocal
                    files.eigenvectorFileName             = "eigenvectors_k_$(k_string).dat"
                    files.eigenvectorShiftedFileName      = "eigenvectors_shifted_k_$(k_string).dat"
                    files.imag_eigenvectorFileName        = "imag_eigenvectors_k_$(k_string).dat"
                    files.imag_eigenvectorShiftedFileName = "imag_eigenvectors_shifted_k_$(k_string).dat"
                    files.frequencyFileName               = "frequencies_k_$(k_string).dat"
                else
                    files.eigenvectorFileName             = "eigenvectors.dat"
                    files.eigenvectorShiftedFileName      = "eigenvectors_shifted.dat"
                    files.frequencyFileName               = "frequencies.dat"
                end

                ########################################
                #                                      #
                # shift potential back to input values #
                #                                      #
                ########################################

                potential.potential = potential.potential .+ potential.shift

                #######################################
                #                                     #
                # convert k-values back to input unit #
                #                                     #
                #######################################

                k = ustrip.(uconvert.(potential.coordsUnit^(-1), k ./ potential.internalElemCoords))

                ############################################################
                #                                                          #
                # print eigenvalues, eigenvectors and frequencies to files #
                #                                                          #
                ############################################################

                printEigenvalues(potential, output, k)
                printEigenvectors(potential, system, output, files, k)
                printFrequencies(potential, system, output, files, k)

                println(i, "/", length(potential.kpoints), " Done")
            end
        end

        ################################################################
        #                                                              #
        # if bandstructure is requested than write band structure file #
        #                                                              #
        ################################################################

        potential.bandStructure && printBandStructure(potential, potential.kpoints)

    end

    ################################
    #                              #
    # print timings of calculation #
    #                              #
    ################################

    files.timingsFile = open(files.timingsFileName, "w")
    show(files.to)
    println()
    show(files.timingsFile, files.to)

    ############
    #          #
    # clean up #
    #          #
    ############

    close(files.timingsFile)
    close(files.logFile)

end