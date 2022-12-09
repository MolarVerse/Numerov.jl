function setupSystem(potential::Potential, system::System1D)

    stencil  = system.stencil
    stencilΔ = system.stencilΔ
    stencil∇ = system.stencil∇
    periodic = system.periodic
    solver   = system.solver
    
    if potential.dimension == 1

        system.n_datapoints = [length(potential.potential)]
    
    elseif potential.dimension == 2
        
        system = System2D()
        
    elseif potential.dimension == 3

        system = System3D()
    end

    system.stencil      = stencil
    system.stencilΔ     = stencilΔ
    system.stencil∇     = stencil∇
    system.periodic     = periodic
    system.n_datapoints = potential.n_datapoints
    system.solver       = solver

    potential.n_kpoints != -1 ? system.reciprocal = true : system.reciprocal = false
    potential.reciprocal = system.reciprocal

    if length(system.periodic) == 1 && potential.dimension != 1
        system.periodic = repeat(system.periodic, potential.dimension)
    end

    !any(system.periodic) && system.reciprocal && (@error "You have defined a number of k-points - this option is only valid in combination with \"periodic = true\""; exit())
    any(system.n_datapoints .< system.stencil) && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    if system.stencil∇ == 0
        system.stencil∇ = system.stencil
    end

    if system.stencilΔ == 0
        system.stencilΔ = system.stencilΔ
    end

    return system
    
end