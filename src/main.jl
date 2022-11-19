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
    system = setupSystem(potential, system)
    buildLaplace(system)
    buildNabla(system)

    checkInput(output)

    isfile("eigenvalues.dat") && rm("eigenvalues.dat")

    if potential.dimension == 1
        for k in potential.kpoints[1]
            solve(potential, system, output, k)

            printEigenvalues(potential, output, k)
            printEigenvectors(potential, system, output, k)
            printFrequencies(potential, system, output, k)
        end
    elseif potential.dimension == 2
        for kx in potential.kpoints[1]
            for ky in potential.kpoints[2]
                solve(potential, system, output, (kx,ky))

                printEigenvalues(potential, output, (kx,ky))
                printEigenvectors(potential, system, output, (kx,ky))
                printFrequencies(potential, system, output, (kx,ky))
            end
        end
    end

end