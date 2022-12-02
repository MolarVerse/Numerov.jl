using TimerOutputs

function numerov(inputFileName::String)

    to = TimerOutput()

    @timeit to "main" begin

    potential = Potential()
    system    = System1D()     #default setup but gets overridden later!
    output    = Output()

    readInputFile(inputFileName)
    checkInput(potential)
    readPotential(potential)

    checkInput(system)
    system = setupSystem(potential, system)
    buildLaplace(system)
    buildNabla(system)

    checkInput(output)

    isfile("eigenvalues.dat") && rm("eigenvalues.dat")

    @timeit to "loop" begin
    for (i, k) in enumerate(potential.kpoints)

        @timeit to "solve" solve(potential, system, output, k, to)

        #shift potential back for output
        potential.potential = potential.potential .+ potential.shift

        @timeit to "print1" printEigenvalues(potential, output, k)
        @timeit to "print2" printEigenvectors(potential, system, output, k)
        @timeit to "print3" printFrequencies(potential, system, output, k)

        println(i, "/", length(potential.kpoints), " Done")
    end
    end

    potential.bandStructure && printBandStructure(potential, potential.kpoints)

    end

    show(to)

end