function buildLaplace(system::System1D)

    n_datapoints = system.n_datapoints[1]
    n_datapoints < system.stencil && (@error "The number of datapoints has at least to be equal to the stencil size!"; exit())

    stencil = zeros(system.stencil)

    if system.stencil == 3

        stencil = [1, -2, 1]

    elseif system.stencil == 5
        
        stencil = [-1/12, 4/3, -5/2, 4/3, -1/12]

    elseif system.stencil == 7

        stencil = [1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90]
        
    elseif system.stencil == 9
        
        stencil = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560]

    elseif system.stencil == 11

        stencil  = [8, -125, 1000, -6000, 42000, -73766, 42000, -6000, 1000, -125, 8]
        stencil /= 25200

    elseif system.stencil == 13

        stencil  = [-50, 864, -7425, 44000, -222750, 1425600, -2480478, 1425600, -222750, 44000, -7425, 864, -50]
        stencil /= 831600

    end

    system.laplace = build_1d_stencil(system, n_datapoints, stencil)

end