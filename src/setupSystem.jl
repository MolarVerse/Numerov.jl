function setupSystem(potential::Potential, system::System1D)

    #######################################################
    #                                                     #
    # setup subtype of system depending on dimensionality #
    #                                                     #
    #######################################################

    stencil  = system.stencil
    stencilΔ = system.stencilΔ
    stencil∇ = system.stencil∇
    solver   = system.solver
    
    if potential.dimension == 1
        system = System1D()
    elseif potential.dimension == 2
        system = System2D()
    elseif potential.dimension == 3
        system = System3D()
    end

    system.stencil      = stencil
    system.stencilΔ     = stencilΔ
    system.stencil∇     = stencil∇
    system.solver       = solver

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

    return system
    
end