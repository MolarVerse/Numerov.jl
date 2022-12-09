function numerov(inputFileName::String)

    potential = Potential()
    system    = System1D()     #default setup but gets overridden later!
    output    = Output()
    files     = Files()

    files.to = TimerOutput()

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

                #shift potential back for output
                potential.potential = potential.potential .+ potential.shift

                printEigenvalues(potential, output, k)
                printEigenvectors(potential, system, output, k)
                printFrequencies(potential, system, output, k)

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