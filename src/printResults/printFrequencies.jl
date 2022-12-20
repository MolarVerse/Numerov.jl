function printFrequencies(potential::Potential, system::System, output::Output, files::Files, k)

    ####################################################
    #                                                  #
    # calculate frequencies as lower triangular matrix #
    #                                                  #
    ####################################################

    output.frequencies = zeros(output.n_eigenvalues-1, output.n_eigenvalues-1)

    for i in 1:output.n_eigenvalues-1
        for j in i+1:output.n_eigenvalues
            output.frequencies[j-1,i] = output.eigenvalues[j] - output.eigenvalues[i]
        end
    end

    ################################################
    #                                              #
    # print frequencies as lower triangular matrix #
    #                                              #
    ################################################

    file = open(files.frequencyFileName, "w")

    @printf(file, "#Lower triangular matrix for frequencies in cm^-1\n")
    
    for i in 1:output.n_eigenvalues-1
        for j in 1:i
            @printf(file, "%12.10lf ", ustrip(uconvert(u"cm^-1", output.frequencies[i,j]*potential.internalElemEnergy / h / c_0)))
        end
        @printf(file, "\n")
    end

    close(file)

end