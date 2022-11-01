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

    checkInput(output)

    solve(potential, system, output)

    printEigenvalues(potential, output)
    printEigenvectors(potential, system, output)
    printFrequencies(potential, output)

end