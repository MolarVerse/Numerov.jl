function printEigenvalues(potential::Potential, output::Output)

    file = open("eigenvalues.dat", "w")

    println(file, "#Eigenvalues given in chosen input unit - $(potential.potentialUnit)")

    k = 0.0 # here k chose as zero but important for band structure
    @printf(file, "%lf ", k)
    for λ in output.eigenvalues
        @printf(file, "%8.6lf ", ustrip(uconvert(potential.potentialUnit, λ*potential.internalElemEnergy)))
    end
    @printf(file, "\n")

    close(file)

end

function printEigenvectors(potential::Potential, system::System, output::Output)
    
    file         = open("eigenvectors.dat", "w")
    file_shifted = open("eigenvectors_shifted.dat", "w")
    
    for i in 1:system.n_datapoints
        for coord in potential.coords
            @printf(file        , "%8.6lf ", ustrip(uconvert(potential.coordsUnit, coord[i]*potential.internalElemCoords)))
            @printf(file_shifted, "%8.6lf ", ustrip(uconvert(potential.coordsUnit, coord[i]*potential.internalElemCoords)))
        end

        @printf(file        , "%8.6lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        @printf(file_shifted, "%8.6lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))

        for (j, ev) in enumerate(output.eigenvectors)
            @printf(file        , "%8.6lf ", ev[i])
            @printf(file_shifted, "%8.6lf ", ev[i] + ustrip(uconvert(potential.potentialUnit, output.eigenvalues[j]*potential.internalElemEnergy)))
        end

        @printf(file        , "\n")
        @printf(file_shifted, "\n")
    end

    close(file)
    close(file_shifted)
    
end

function printFrequencies(potential::Potential, output::Output)

    for i in 1:output.n_eigenvalues-1
        for j in i+1:output.n_eigenvalues
            output.frequencies[j-1,i] = output.eigenvalues[j] - output.eigenvalues[i]
        end
    end

    file = open("frequencies.dat", "w")

    @printf(file, "#Lower triangular matrix for frequencies in cm^-1\n")
    
    for i in 1:output.n_eigenvalues-1
        for j in 1:i
            @printf(file, "%8.6lf ", ustrip(uconvert(u"cm^-1", output.frequencies[i,j]*unit(potential.internalElemEnergy) / h / c_0)))
        end
        @printf(file, "\n")
    end

    close(file)

end