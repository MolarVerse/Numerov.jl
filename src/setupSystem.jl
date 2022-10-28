function setupSystem(potential::Potential, system::System1D)

    system.n_datapoints = length(potential.potential)
    
end