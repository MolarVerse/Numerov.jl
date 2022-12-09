function numerov(inputFileName::String)

    potential = Potential()
    system    = System1D()     #default setup but gets overridden later!
    output    = Output()
    files     = Files()

    files.to = TimerOutput()

    @timeit files.to "main" begin

        readInputFile(inputFileName)
        
        checkInput(potential)
        checkInput(files)
        checkInput(output)

        readPotential(potential)

        checkInput(system)
        system = setupSystem(potential, system)
        buildΔ(system)
        build∇(system)


        isfile("eigenvalues.dat") && rm("eigenvalues.dat")

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
    show(files.timingsFile, files.to)

    close(files.timingsFile)
    close(files.logFile)

end