function printEigenvalues(potential::Potential, output::Output, k)

    file = open("eigenvalues.dat", "a")

    filesize("eigenvalues.dat") == 0 && println(file, "#Eigenvalues given in chosen input unit - $(potential.potentialUnit)")

    if potential.reciprocal
        for i in 1:potential.dimension
            @printf(file, "%lf ", k[i])
        end
    end
    for λ in output.eigenvalues
        @printf(file, "%20.14lf ", ustrip(uconvert(potential.potentialUnit, λ*potential.internalElemEnergy)))
    end
    @printf(file, "\n")

    close(file)

end

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

        @printf(file        , "%8.6lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        @printf(file_shifted, "%8.6lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        if system.reciprocal
            @printf(imag_file        , "%8.6lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
            @printf(imag_file_shifted, "%8.6lf ", ustrip(uconvert(potential.potentialUnit, potential.potential[i]*potential.internalElemEnergy)))
        end

        for (j, ev) in enumerate(output.eigenvectors)
            @printf(file        , "%8.6lf ", real(ev[i]))
            @printf(file_shifted, "%8.6lf ", real(ev[i]) + ustrip(uconvert(potential.potentialUnit, (output.eigenvalues[j]+ potential.shift)*potential.internalElemEnergy)))
            if system.reciprocal
                @printf(imag_file        , "%8.6lf ", imag(ev[i]))
                @printf(imag_file_shifted, "%8.6lf ", imag(ev[i]) + ustrip(uconvert(potential.potentialUnit, (output.eigenvalues[j]+ potential.shift)*potential.internalElemEnergy)))
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

function printFrequencies(potential::Potential, system::System, output::Output, k)

    for i in 1:output.n_eigenvalues-1
        for j in i+1:output.n_eigenvalues
            output.frequencies[j-1,i] = output.eigenvalues[j] - output.eigenvalues[i]
        end
    end

    if system.reciprocal
        file = open("frequencies_k=$(k).dat", "w")
    else
        file = open("frequencies.dat", "w")
    end

    @printf(file, "#Lower triangular matrix for frequencies in cm^-1\n")
    
    for i in 1:output.n_eigenvalues-1
        for j in 1:i
            @printf(file, "%8.6lf ", ustrip(uconvert(u"cm^-1", output.frequencies[i,j]*potential.internalElemEnergy / h / c_0)))
        end
        @printf(file, "\n")
    end

    close(file)

end

function printBandStructure(potential::Potential, k_points)

    file = open("bandstructure.dat", "w")

    data = readdlm("eigenvalues.dat"; skipstart=1)

    a = copy(k_points)
    b = copy(k_points)

    pushfirst!(a, Tuple(zeros(length(k_points[1]))))
    push!(b, Tuple(zeros(length(k_points[1]))))

    diff = (norm.([b[i] .- a[i] for i in eachindex(a)]))[1:end-1]

    brioullin_path = 0.0

    for (i, spacing) in enumerate(diff)
        brioullin_path += spacing
        @printf(file, "%lf ", brioullin_path)
        for j in length(k_points[1])+1:length(data[1,:])
            @printf(file, "%20.14lf ", data[i,j])
        end
        @printf(file, "\n")
        if i == length(diff) - potential.n_kpoints && potential.dimension == 3
            break
        end
    end

    for (i, spacing) in enumerate(diff)
        if i <= length(diff) - potential.n_kpoints || potential.dimension != 3
            continue
        end
        @printf(file, "%lf ", brioullin_path)
        for j in length(k_points[1])+1:length(data[1,:])
            @printf(file, "%20.14lf ", data[i,j])
        end
        @printf(file, "\n")
        brioullin_path += spacing
    end

    close(file)    
end