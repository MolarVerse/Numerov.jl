function setupSystem(potential::Potential, system::System)

    system.periodic     = potential.periodic
    system.n_datapoints = potential.n_datapoints
    system.reciprocal   = potential.reciprocal

    ##############################################################################################
    #                                                                                            #
    # check if system is periodic if k-points are set and minimal length of potential datapoints #
    #                                                                                            #
    ##############################################################################################

    system.reciprocal && !any(system.periodic) && (@error "You have defined a number of k-points - this option is only valid in combination with \"periodic = true\""; exit())
    any(system.n_datapoints .< system.stencil) && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    ###############################################################
    #                                                             #
    # set stencil for laplace and nabla if not defined seperately #
    #                                                             #
    ###############################################################

    system.stencil∇ == 0 && (system.stencil∇ = system.stencil)
    system.stencilΔ == 0 && (system.stencilΔ = system.stencil)

end