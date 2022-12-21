function printEigenvalues(potential::Potential, output::Output, files::Files, k)

    #########################################################################
    #                                                                       #
    # open file and check if size is zero in order to print info about unit #
    #                                                                       #
    #########################################################################

    file = open(files.eigenvalueFileName, "a")

    filesize(files.eigenvalueFileName) == 0 && println(file, "#Eigenvalues given in chosen input unit - $(potential.potentialUnit)")

    ##########################################################################
    #                                                                        #
    # print k-points if system is requested to calculate in reciprocal space #
    #                                                                        #
    ##########################################################################

    if potential.reciprocal
        for i in 1:potential.dimension
            @printf(file, "%lf ", k[i])
        end
    end

    #################################################
    #                                               #
    # print all requested eigenvalues in input unit #
    #                                               #
    #################################################

    for λ in output.eigenvalues
        @printf(file, "%20.14lf ", ustrip(uconvert(potential.potentialUnit, λ*potential.internalElemEnergy)))
    end
    @printf(file, "\n")

    ############
    #          #
    # clean up #
    #          #
    ############

    close(file)

end