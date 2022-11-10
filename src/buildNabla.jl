function buildNabla(system::System1D)

    n_datapoints = system.n_datapoints
    n_datapoints < system.stencil && (@error "The number of datapoints has at least to be equal to the stencil size!"; exit())

    !system.bandStructure && (system.Δ = sparse(zeros(n_datapoints, n_datapoints)); return)

    if system.stencil == 3

        system.Δ = spdiagm(
            -1 => ones(n_datapoints-1)*(-1/2),
             1 => ones(n_datapoints-1)*( 1/2)
        )

        if system.periodic  #remove this periodic argument

            system.Δ += spdiagm(
                 n_datapoints-1 => ones(1)*(-1/2),
                -n_datapoints+1 => ones(1)*( 1/2)
            )

        end

    elseif system.stencil == 5
        
        system.Δ = spdiagm(
            -2 => ones(n_datapoints-2)*( 1/12),
            -1 => ones(n_datapoints-1)*(-2/3 ),
             1 => ones(n_datapoints-1)*( 2/3 ),
             2 => ones(n_datapoints-2)*(-1/12)
        )

        if system.periodic

            system.Δ += spdiagm(
                 n_datapoints-2 => ones(2)*( 1/12),
                 n_datapoints-1 => ones(1)*(-2/3 ),
                -n_datapoints+1 => ones(1)*( 2/3 ),
                -n_datapoints+2 => ones(2)*(-1/12)
            )
            
        end

    elseif system.stencil == 7
        
        system.Δ = spdiagm(
            -3 => ones(n_datapoints-3)*( -1/60),
            -2 => ones(n_datapoints-2)*(  3/20),
            -1 => ones(n_datapoints-1)*( -3/4 ),
             1 => ones(n_datapoints-1)*(  3/4 ),
             2 => ones(n_datapoints-2)*( -3/20),
             3 => ones(n_datapoints-3)*(  1/60)
        )

        if system.periodic

            system.Δ += spdiagm(
                 n_datapoints-3 => ones(3)*( -1/60),
                 n_datapoints-2 => ones(2)*(  3/20),
                 n_datapoints-1 => ones(1)*( -3/4 ),
                -n_datapoints+1 => ones(1)*(  3/4 ),
                -n_datapoints+2 => ones(2)*( -3/20),
                -n_datapoints+3 => ones(3)*(  1/60)
            )
            
        end

    elseif system.stencil == 9
        
        system.Δ = spdiagm(
            -4 => ones(n_datapoints-4)*(   1/280),
            -3 => ones(n_datapoints-3)*(  -4/105),
            -2 => ones(n_datapoints-2)*(   1/5  ),
            -1 => ones(n_datapoints-1)*(  -4/5  ),
             1 => ones(n_datapoints-1)*(   4/5  ),
             2 => ones(n_datapoints-2)*(  -1/5  ),
             3 => ones(n_datapoints-3)*(   4/105),
             4 => ones(n_datapoints-4)*(  -1/280)
        )

        if system.periodic

            system.Δ += spdiagm(
                 n_datapoints-4 => ones(4)*(   1/280),
                 n_datapoints-3 => ones(3)*(  -4/105),
                 n_datapoints-2 => ones(2)*(   1/5  ),
                 n_datapoints-1 => ones(1)*(  -4/5  ),
                -n_datapoints+1 => ones(1)*(   4/5  ),
                -n_datapoints+2 => ones(2)*(  -1/5  ),
                -n_datapoints+3 => ones(3)*(   4/105),
                -n_datapoints+4 => ones(4)*(  -1/280)
            )
            
        end

    elseif system.stencil == 11
        
        system.Δ = spdiagm(
            -5 => ones(n_datapoints-5)*(   -2/25200),
            -4 => ones(n_datapoints-4)*(   25/25200),
            -3 => ones(n_datapoints-3)*( -150/25200),
            -2 => ones(n_datapoints-2)*(  600/25200),
            -1 => ones(n_datapoints-1)*(-2100/25200),
             1 => ones(n_datapoints-1)*( 2100/25200),
             2 => ones(n_datapoints-2)*( -600/25200),
             3 => ones(n_datapoints-3)*(  150/25200),
             4 => ones(n_datapoints-4)*(  -25/25200),
             5 => ones(n_datapoints-5)*(    2/25200)
        )

        if system.periodic

            system.Δ += spdiagm(
                 n_datapoints-5 => ones(5)*(   -2/25200),
                 n_datapoints-4 => ones(4)*(   25/25200),
                 n_datapoints-3 => ones(3)*( -150/25200),
                 n_datapoints-2 => ones(2)*(  600/25200),
                 n_datapoints-1 => ones(1)*(-2100/25200),
                -n_datapoints+1 => ones(1)*( 2100/25200),
                -n_datapoints+2 => ones(2)*( -600/25200),
                -n_datapoints+3 => ones(3)*(  150/25200),
                -n_datapoints+4 => ones(4)*(  -25/25200),
                -n_datapoints+5 => ones(5)*(    2/25200)
            )
            
        end

    elseif system.stencil == 13

        system.bandStructure && (@error "13-point stencil is not yet implemented for bandstructure calculation!"; exit())

    end

end