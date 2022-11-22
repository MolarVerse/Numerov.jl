function buildLaplace(system::System1D)

    n_datapoints = system.n_datapoints
    n_datapoints < system.stencil && (@error "The number of datapoints has at least to be equal to the stencil size!"; exit())

    if system.stencil == 3

        system.laplace = spdiagm(
            -1 => ones(n_datapoints-1),
             0 => ones(n_datapoints  )*(-2),
             1 => ones(n_datapoints-1)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-1 => ones(1),
                -n_datapoints+1 => ones(1)
            )

        end

    elseif system.stencil == 5
        
        system.laplace = spdiagm(
            -2 => ones(n_datapoints-2)*(-1/12),
            -1 => ones(n_datapoints-1)*( 4/3 ),
             0 => ones(n_datapoints  )*(-5/2 ),
             1 => ones(n_datapoints-1)*( 4/3 ),
             2 => ones(n_datapoints-2)*(-1/12)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-2 => ones(2)*(-1/12),
                 n_datapoints-1 => ones(1)*( 4/3),
                -n_datapoints+1 => ones(1)*( 4/3),
                -n_datapoints+2 => ones(2)*(-1/12)
            )
            
        end

    elseif system.stencil == 7
        
        system.laplace = spdiagm(
            -3 => ones(n_datapoints-3)*(  1/90),
            -2 => ones(n_datapoints-2)*( -3/20),
            -1 => ones(n_datapoints-1)*(  3/2 ),
             0 => ones(n_datapoints  )*(-49/18),
             1 => ones(n_datapoints-1)*(  3/2 ),
             2 => ones(n_datapoints-2)*( -3/20),
             3 => ones(n_datapoints-3)*(  1/90)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-3 => ones(3)*(  1/90),
                 n_datapoints-2 => ones(2)*( -3/20),
                 n_datapoints-1 => ones(1)*(  3/2 ),
                -n_datapoints+1 => ones(1)*(  3/2 ),
                -n_datapoints+2 => ones(2)*( -3/20),
                -n_datapoints+3 => ones(3)*(  1/90)
            )
            
        end

    elseif system.stencil == 9
        
        system.laplace = spdiagm(
            -4 => ones(n_datapoints-4)*(  -1/560),
            -3 => ones(n_datapoints-3)*(   8/315),
            -2 => ones(n_datapoints-2)*(  -1/5  ),
            -1 => ones(n_datapoints-1)*(   8/5  ),
             0 => ones(n_datapoints  )*(-205/72 ),
             1 => ones(n_datapoints-1)*(   8/5  ),
             2 => ones(n_datapoints-2)*(  -1/5  ),
             3 => ones(n_datapoints-3)*(   8/315),
             4 => ones(n_datapoints-4)*(  -1/560)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-4 => ones(4)*(  -1/560),
                 n_datapoints-3 => ones(3)*(   8/315),
                 n_datapoints-2 => ones(2)*(  -1/5  ),
                 n_datapoints-1 => ones(1)*(   8/5  ),
                -n_datapoints+1 => ones(1)*(   8/5  ),
                -n_datapoints+2 => ones(2)*(  -1/5  ),
                -n_datapoints+3 => ones(3)*(   8/315),
                -n_datapoints+4 => ones(4)*(  -1/560)
            )
            
        end

    elseif system.stencil == 11
        
        system.laplace = spdiagm(
            -5 => ones(n_datapoints-5)*(     8/25200),
            -4 => ones(n_datapoints-4)*(  -125/25200),
            -3 => ones(n_datapoints-3)*(  1000/25200),
            -2 => ones(n_datapoints-2)*( -6000/25200),
            -1 => ones(n_datapoints-1)*( 42000/25200),
             0 => ones(n_datapoints  )*(-73766/25200),
             1 => ones(n_datapoints-1)*( 42000/25200),
             2 => ones(n_datapoints-2)*( -6000/25200),
             3 => ones(n_datapoints-3)*(  1000/25200),
             4 => ones(n_datapoints-4)*(  -125/25200),
             5 => ones(n_datapoints-5)*(     8/25200)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-5 => ones(5)*(     8/25200),
                 n_datapoints-4 => ones(4)*(  -125/25200),
                 n_datapoints-3 => ones(3)*(  1000/25200),
                 n_datapoints-2 => ones(2)*( -6000/25200),
                 n_datapoints-1 => ones(1)*( 42000/25200),
                -n_datapoints+1 => ones(1)*( 42000/25200),
                -n_datapoints+2 => ones(2)*( -6000/25200),
                -n_datapoints+3 => ones(3)*(  1000/25200),
                -n_datapoints+4 => ones(4)*(  -125/25200),
                -n_datapoints+5 => ones(5)*(     8/25200)
            )
            
        end

    elseif system.stencil == 13

        system.laplace = spdiagm(
            -6 => ones(n_datapoints-6)*(     -50/831600),
            -5 => ones(n_datapoints-5)*(     864/831600),
            -4 => ones(n_datapoints-4)*(   -7425/831600),
            -3 => ones(n_datapoints-3)*(   44000/831600),
            -2 => ones(n_datapoints-2)*( -222750/831600),
            -1 => ones(n_datapoints-1)*( 1425600/831600),
             0 => ones(n_datapoints  )*(-2480478/831600),
             1 => ones(n_datapoints-1)*( 1425600/831600),
             2 => ones(n_datapoints-2)*( -222750/831600),
             3 => ones(n_datapoints-3)*(   44000/831600),
             4 => ones(n_datapoints-4)*(   -7425/831600),
             5 => ones(n_datapoints-5)*(     864/831600),
             6 => ones(n_datapoints-6)*(     -50/831600)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-6 => ones(6)*(     -50/831600),
                 n_datapoints-5 => ones(5)*(     864/831600),
                 n_datapoints-4 => ones(4)*(   -7425/831600),
                 n_datapoints-3 => ones(3)*(   44000/831600),
                 n_datapoints-2 => ones(2)*( -222750/831600),
                 n_datapoints-1 => ones(1)*( 1425600/831600),
                -n_datapoints+1 => ones(1)*( 1425600/831600),
                -n_datapoints+2 => ones(2)*( -222750/831600),
                -n_datapoints+3 => ones(3)*(   44000/831600),
                -n_datapoints+4 => ones(4)*(   -7425/831600),
                -n_datapoints+5 => ones(5)*(     864/831600),
                -n_datapoints+6 => ones(6)*(     -50/831600)
            )
            
        end

    end

end