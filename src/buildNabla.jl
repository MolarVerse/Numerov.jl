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

function buildNabla(system::System2D)

    n_datapoints = system.n_datapoints

    !system.bandStructure && (system.Δ = sparse(zeros(prod(n_datapoints), prod(n_datapoints))); return)
    
    n_datapoints[1] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())
    n_datapoints[2] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    system.Δ   = zeros(n_datapoints[1]*n_datapoints[2], n_datapoints[1]*n_datapoints[2])
    stencil    = zeros(system.stencil, system.stencil)
    stencil_1d = zeros(system.stencil)

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

    stencil[:,system.stencil÷2+1] = stencil_1d
    stencil[system.stencil÷2+1,:] = stencil_1d

    for x_1 in 1:n_datapoints[1]
        for stencil_index in 1:system.stencil
    
            super_matrix_index = x_1 - system.stencil ÷ 2 - 1 + stencil_index
    
            if super_matrix_index < 1 && system.periodic == false
                continue
            elseif super_matrix_index < 1
                super_matrix_index += n_datapoints[1]
            end
    
            matrix = zeros(n_datapoints[2], n_datapoints[2])
    
            for j in 1:system.stencil
    
                sub_matrix_index = j - system.stencil ÷ 2 - 1
    
                if sub_matrix_index == 0
                    matrix += spdiagm(sub_matrix_index => ones(n_datapoints[2] - sub_matrix_index) * stencil[stencil_index, j])
                else
                    matrix += spdiagm(sub_matrix_index => ones(n_datapoints[2] - abs(sub_matrix_index)) * stencil[stencil_index, j])
                end
    
                if system.periodic && sub_matrix_index != 0
                    matrix += spdiagm(sign(sub_matrix_index)*n_datapoints[2] - sub_matrix_index => ones(abs(sub_matrix_index)) * sign(-sub_matrix_index)*stencil[stencil_index, j])
                end
            end
    
            i_1 = 1 + (x_1 - 1) * n_datapoints[2]
            j_1 = x_1 * n_datapoints[2]
    
            i_2 = 1 + (super_matrix_index - 1) * n_datapoints[2]
            j_2 = (super_matrix_index) * n_datapoints[2]
    
            if(i_2 > n_datapoints[2]*n_datapoints[1])
                continue
            end
    
            system.Δ[i_1:j_1, i_2:j_2] += matrix #add system.laplace
    
        end
    end


end