function buildNabla(system::System1D) #combine these two functions!

    n_datapoints = system.n_datapoints
    n_datapoints < system.stencil && (@error "The number of datapoints has at least to be equal to the stencil size!"; exit())

    !system.bandStructure && (system.Δ = sparse(zeros(n_datapoints, n_datapoints)); return)

    stencil = get_1d_stencil(system)

    system.Δ = build_1d_stencil(system, n_datapoints, stencil)

end

function buildNabla(system::System2D)

    n_datapoints = system.n_datapoints

    !system.bandStructure && (system.Δ = sparse(zeros(prod(n_datapoints), prod(n_datapoints))); return)
    
    n_datapoints[1] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())
    n_datapoints[2] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    system.Δ   = zeros(n_datapoints[1]*n_datapoints[2], n_datapoints[1]*n_datapoints[2])
    stencil    = zeros(system.stencil, system.stencil)
    stencil_1d = get_1d_stencil(system)

    stencil[:,system.stencil÷2+1] = stencil_1d
    stencil[system.stencil÷2+1,:] = stencil_1d

    system.Δ = build_2d_stencil(system, n_datapoints, stencil)

end

function get_1d_stencil(system)

    stencil = zeros(system.stencil)

    if system.stencil == 3

        stencil_1d  = [-1, 0, 1]
        stencil_1d /= 2.0

    elseif system.stencil == 5
        
        stencil_1d  = [1, -8, 0, 8, -1]
        stencil_1d /= 12.0

    elseif system.stencil == 7
        
        stencil_1d  = [0, -9, -45, 0, 45, -9, 1]
        stencil_1d /= 60.0

    elseif system.stencil == 9
        
        stencil_1d  = [3, -32, 168, -672, 0, 672, -168, 32, -3]
        stencil_1d /= 840.0

    elseif system.stencil == 11
        
        stencil_1d  = [-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2]
        stencil_1d /= 25200.0

    elseif system.stencil == 13

        system.bandStructure && (@error "13-point stencil is not yet implemented for bandstructure calculation!"; exit())

    end   
    
    return stencil
end