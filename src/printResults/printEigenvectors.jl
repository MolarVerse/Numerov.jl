function printEigenvectors(potential::Potential, system::System, output::Output, k)   ###### probably not normalized!!!!!!!!!!!

    if system.reciprocal
        file              = open("eigenvectors_k=$(k).dat", "w")
        file_shifted      = open("eigenvectors_shifted_k=$(k).dat", "w")
        imag_file         = open("imag_eigenvectors_k=$(k).dat", "w")
        imag_file_shifted = open("imag_eigenvectors_shifted_k=$(k).dat", "w")
    else
        file         = open("eigenvectors.dat", "w")
        file_shifted = open("eigenvectors_shifted.dat", "w")
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