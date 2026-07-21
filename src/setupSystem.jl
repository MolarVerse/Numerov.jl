function setupSystem(potential::Potential, system::System)

    system.periodic     = potential.periodic
    system.n_datapoints = potential.n_datapoints
    system.reciprocal   = potential.reciprocal

    ##############################################################################################
    #                                                                                            #
    # check if system is periodic if k-points are set and minimal length of potential datapoints #
    #                                                                                            #
    ##############################################################################################

    system.reciprocal && !any(system.periodic) && throw(ArgumentError("You have defined a number of k-points - this option is only valid in combination with \"periodic = true\""))
    any(system.n_datapoints .< system.stencil) && throw(ArgumentError("The number of datapoints in each dimension has at least to be equal to the stencil size!"))

    # solveWrapper only ever builds a real Hamiltonian for lobpcg, so a
    # reciprocal (k-point) run has to be rejected here - before any output
    # file is written or an existing eigenvalues.dat is removed - rather than
    # deep inside solve() on the first k-point
    system.reciprocal && system.solver == LOBPCG && throw(ArgumentError("the lobpcg solver supports non-periodic (real symmetric) problems - use arpack or krylov for periodic k-point runs"))

    ###############################################################
    #                                                             #
    # set stencil for laplace and nabla if not defined seperately #
    #                                                             #
    ###############################################################

    system.stencil∇ == 0 && (system.stencil∇ = system.stencil)
    system.stencilΔ == 0 && (system.stencilΔ = system.stencil)

end