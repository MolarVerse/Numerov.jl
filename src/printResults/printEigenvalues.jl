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