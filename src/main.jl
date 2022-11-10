function numerov(inputFileName::String)

    potential = Potential()
    system    = System1D()     #default setup but gets overridden later!
    output    = Output()

    readInputFile(inputFileName)
    checkInput(potential)
    readPotential(potential)

    # if potential.dimension == 1
    #     potential.periodic ? exit() : (system = System1D())
    # end

    checkInput(system)
    setupSystem(potential, system)
    buildLaplace(system)
    buildNabla(system)

    checkInput(output)

    isfile("eigenvalues.dat") && rm("eigenvalues.dat")

    for k in potential.kpoints
        solve(potential, system, output, k)

        printEigenvalues(potential, output, k)
        printEigenvectors(potential, system, output, k)
        printFrequencies(potential, system, output, k)
    end

end