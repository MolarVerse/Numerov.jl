function numerov(inputFileName::String)

    potential = Potential()
    system    = System1D()     #default setup but gets overridden later!
    output    = Output()
    files     = Files()

    files.to = TimerOutput()

    [inputDictionary[key] = "" for key in keys(inputDictionary)] #to reset dict to default values if calculation started in same repl session

    @timeit files.to "main" begin

        readInputFile(inputFileName)
        
        checkInput(potential)
        checkInput(system)
        checkInput(files)
        checkInput(output)

        readPotential(potential)

        system = setupSystem(potential, system)

        buildΔ(system)
        build∇(system)

        isfile("eigenvalues.dat") && rm("eigenvalues.dat")

        println(files.logFile, "Non zeros = ", length(system.Δ.nzval))
        println(files.logFile, "zeros     = ", length(system.Δ) - length(system.Δ.nzval))
        println(files.logFile, "Sparsity  = ", length(system.Δ.nzval) / length(system.Δ))

        @timeit files.to "loop" begin
            for (i, k) in enumerate(potential.kpoints)

                @timeit files.to "solve" solve(potential, system, output, k, files)

                k = k .* sqrt.(potential.mass)

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

                #shift potential back for output
                potential.potential = potential.potential .+ potential.shift

                k = ustrip.(uconvert.(potential.coordsUnit^(-1), k ./ potential.internalElemCoords))

                printEigenvalues(potential, output, k)
                printEigenvectors(potential, system, output, files, k)
                printFrequencies(potential, system, output, files, k)

                println(i, "/", length(potential.kpoints), " Done")
            end
        end

        potential.bandStructure && printBandStructure(potential, potential.kpoints)

    end

    files.timingsFile = open(files.timingsFileName, "w")
    show(files.to)
    println()
    show(files.timingsFile, files.to)

    close(files.timingsFile)
    close(files.logFile)

end