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
            solve(potential, system, output, k, to)

            printEigenvalues(potential, output, k)
            printEigenvectors(potential, system, output, k)
            printFrequencies(potential, system, output, k)
        end
    elseif potential.dimension == 2

        @timeit to "loop" begin
            for (i, kx) in enumerate(potential.kpoints[1])
                for (j, ky) in enumerate(potential.kpoints[2])
    
                    @timeit to "solve" solve(potential, system, output, (kx,ky), to)
    
                    @timeit to "print1" printEigenvalues(potential, output, (kx,ky))
                    @timeit to "print2" printEigenvectors(potential, system, output, (kx,ky))
                    @timeit to "print3" printFrequencies(potential, system, output, (kx,ky))
    
                    println((i-1)*length(potential.kpoints[1]) + j, "/", potential.n_kpoints*potential.n_kpoints, " Done")
    
                end
            end
        end
    elseif potential.dimension == 3
        @timeit to "loop" begin
            for (i, kx) in enumerate(potential.kpoints[1])
                for (j, ky) in enumerate(potential.kpoints[2])
                    for (j, kz) in enumerate(potential.kpoints[3])
    
                    @timeit to "solve" solve(potential, system, output, (kx,ky,kz), to)
    
                    @timeit to "print1" printEigenvalues(potential, output, (kx,ky, kz))
                    @timeit to "print2" printEigenvectors(potential, system, output, (kx,ky, kz))
                    @timeit to "print3" printFrequencies(potential, system, output, (kx,ky, kz))
       
                    end

                    println((i-1)*length(potential.kpoints[1]) + j, "/", potential.n_kpoints*potential.n_kpoints, " Done")

                end
            end
        end
    end

    end

    show(to)

end