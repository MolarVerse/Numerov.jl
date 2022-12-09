function buildNabla(system::System1D) #combine these two functions!

    n_datapoints = system.n_datapoints[1]
    n_datapoints < system.stencil && (@error "The number of datapoints has at least to be equal to the stencil size!"; exit())

    !system.reciprocal && (system.∇ = spzeros(n_datapoints, n_datapoints); return)

    stencil = get_1d_stencil(system)

    system.∇ = build_1d_stencil(system, n_datapoints, stencil)

end

function buildNabla(system::System2D)

    n_datapoints = system.n_datapoints

    !system.reciprocal && (system.∇ = spzeros(prod(n_datapoints), prod(n_datapoints)); return)
    
    n_datapoints[1] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())
    n_datapoints[2] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    system.∇   = spzeros(prod(n_datapoints), prod(n_datapoints))
    stencil    = zeros(system.stencil, system.stencil)
    stencil_1d = get_1d_stencil(system)

    stencil[:,system.stencil÷2+1] = stencil_1d
    stencil[system.stencil÷2+1,:] = stencil_1d

    system.∇ = build_2d_stencil(system, n_datapoints, stencil)
end

function buildNabla(system::System3D)

    n_datapoints = system.n_datapoints

    !system.reciprocal && (system.∇ = spzeros(prod(n_datapoints), prod(n_datapoints)); return)
    
    n_datapoints[1] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())
    n_datapoints[2] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())
    n_datapoints[3] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    system.∇   = spzeros(prod(n_datapoints), prod(n_datapoints))
    stencil    = zeros(system.stencil, system.stencil, system.stencil)
    stencil_1d = get_1d_stencil(system)
# 
    stencil[:                 ,system.stencil÷2+1,system.stencil÷2+1] = stencil_1d
    stencil[system.stencil÷2+1,:                 ,system.stencil÷2+1] = stencil_1d
    stencil[system.stencil÷2+1,system.stencil÷2+1,:                 ] = stencil_1d
# 
    system.∇ = build_3d_stencil(system, n_datapoints, stencil)

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
        stencil_1d /= 2520.0

    elseif system.stencil == 13

        system.reciprocal && (@error "13-point stencil is not yet implemented for reciprocal calculation!"; exit())

    end   
    
    return stencil_1d
end