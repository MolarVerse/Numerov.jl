function setupSystem(potential::Potential, system::System1D)

    if potential.dimension == 1

        system.n_datapoints = length(potential.potential)
    
    elseif potential.dimension == 2
        
        stencil = system.stencil
        periodic = system.periodic
        
        system = System2D()

        system.stencil = stencil
        system.periodic = periodic
        system.n_datapoints = potential.n_datapoints
    elseif potential.dimension == 3

        stencil = system.stencil
        periodic = system.periodic

        system = System3D()

        system.stencil = stencil
        system.periodic = periodic
        system.n_datapoints = potential.n_datapoints
    end

    system.Δ = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))

    potential.n_kpoints != -1 ? system.reciprocal = true : system.reciprocal = false
    potential.reciprocal = system.reciprocal

    !system.periodic && system.reciprocal && (@error "You have defined a number of k-points - this option is only valid in combination with \"periodic = true\""; exit())

    return system
    
end