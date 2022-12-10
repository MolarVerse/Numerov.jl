function printEigenvectors(potential::Potential, system::System, output::Output, files::Files, k)

    file = open(files.eigenvectorFileName, "w")
    file_shifted = open(files.eigenvectorShiftedFileName, "w")

    if system.reciprocal
        imag_file = open(files.imag_eigenvectorFileName, "w")
        imag_file_shifted = open(files.imag_eigenvectorShiftedFileName, "w")
    end

    #think of a clever way to handle shifted input potential for output
    
    for i in 1:prod(system.n_datapoints)
        for coord in potential.coords
            @printf(file        , "%8.6lf ", ustrip(uconvert(potential.coordsUnit, coord[i]*potential.internalElemCoords)))
            @printf(file_shifted, "%8.6lf ", ustrip(uconvert(potential.coordsUnit, coord[i]*potential.internalElemCoords)))
            if system.reciprocal
                @printf(imag_file        , "%8.6lf ", ustrip(uconvert(potential.coordsUnit, coord[i]*potential.internalElemCoords)))
                @printf(imag_file_shifted, "%8.6lf ", ustrip(uconvert(potential.coordsUnit, coord[i]*potential.internalElemCoords)))
            end
        end

        @printf(file        , "%12.10lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        @printf(file_shifted, "%12.10lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        if system.reciprocal
            @printf(imag_file        , "%12.10lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
            @printf(imag_file_shifted, "%12.10lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        end

        for (j, ev) in enumerate(output.eigenvectors)
            @printf(file        , "%12.10lf ", real(ev[i]))
            @printf(file_shifted, "%12.10lf ", real(ev[i]) + ustrip(uconvert(potential.potentialUnit, (output.eigenvalues[j]+ potential.shift)*potential.internalElemEnergy)))
            if system.reciprocal
                @printf(imag_file        , "%12.10lf ", imag(ev[i]))
                @printf(imag_file_shifted, "%12.10lf ", imag(ev[i]) + ustrip(uconvert(potential.potentialUnit, (output.eigenvalues[j]+ potential.shift)*potential.internalElemEnergy)))
            end
        end

        @printf(file        , "\n")
        @printf(file_shifted, "\n")
        if system.reciprocal
            @printf(imag_file        , "\n")
            @printf(imag_file_shifted, "\n")
        end
    end

    close(file)
    close(file_shifted)
    if system.reciprocal
        close(imag_file)
        close(imag_file_shifted)
    end
end