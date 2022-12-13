function build∇(system::System1D) #combine these two functions!

    stencil  = get_1d_stencil(system)
    system.∇ = build_1d_stencil(system, stencil, system.stencil∇)

end

function build∇(system::System2D)

    system.∇   = spzeros(prod(system.n_datapoints), prod(system.n_datapoints))
    stencil    = zeros(system.stencil∇, system.stencil∇)
    stencil_1d = get_1d_stencil(system)

    stencil[:,system.stencil∇÷2+1] = stencil_1d
    stencil[system.stencil∇÷2+1,:] = stencil_1d

    system.∇ = build_2d_stencil(system, system.n_datapoints, stencil, system.stencil∇)
end

function build∇(system::System3D)

    system.∇   = spzeros(prod(system.n_datapoints), prod(system.n_datapoints))
    stencil    = zeros(system.stencil∇, system.stencil∇, system.stencil∇)
    stencil_1d = get_1d_stencil(system)
 
    stencil[:                 ,system.stencil∇÷2+1,system.stencil∇÷2+1] = stencil_1d
    stencil[system.stencil∇÷2+1,:                 ,system.stencil∇÷2+1] = stencil_1d
    stencil[system.stencil∇÷2+1,system.stencil∇÷2+1,:                 ] = stencil_1d
 
    system.∇ = build_3d_stencil(system, system.n_datapoints, stencil, system.stencil∇)

end

function get_1d_stencil(system)

    if system.stencil∇ == 3

        stencil_1d  = [-1, 0, 1] ./ 2.0

    elseif system.stencil∇ == 5
        
        stencil_1d  = [1, -8, 0, 8, -1] ./ 12.0

    elseif system.stencil∇ == 7
        
        stencil_1d  = [0, -9, -45, 0, 45, -9, 1] ./ 60.0

    elseif system.stencil∇ == 9
        
        stencil_1d  = [3, -32, 168, -672, 0, 672, -168, 32, -3] ./ 840.0

    elseif system.stencil∇ == 11
        
        stencil_1d  = [-2, 25, -150, 600, -2100, 0, 2100, -600, 150, -25, 2] ./ 2520.0

    elseif system.stencil∇ == 13

        system.reciprocal && (@error "13-point stencil is not yet implemented for reciprocal calculation!"; exit())

    end   
    
    return stencil_1d
end

function build∇_k(potential::Potential, system::System, k)
    ∇ = spzeros(prod(potential.n_datapoints), prod(potential.n_datapoints))

    if potential.reciprocal
        if potential.dimension == 1
            ∇ = system.∇*k[1]
        elseif potential.dimension == 2
            if k[2] != 0.0
                ∇ = system.∇ * k[2]
                for i in 1:potential.n_datapoints[1]
                    ∇[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] *= k[1]/k[2]
                end
            else
                for i in 1:potential.n_datapoints[1]
                    ∇[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] = system.∇[(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2]] * k[1]
                end
            end

        elseif potential.dimension == 3

            stencil    = zeros(system.stencil, system.stencil, system.stencil)

            stencil[:                 ,system.stencil÷2+1,system.stencil÷2+1] = ones(system.stencil)*k[1]
            stencil[system.stencil÷2+1,:                 ,system.stencil÷2+1] = ones(system.stencil)*k[2]
            stencil[system.stencil÷2+1,system.stencil÷2+1,:                 ] = ones(system.stencil)*k[3]
            stencil[system.stencil÷2+1,system.stencil÷2+1,system.stencil÷2+1] = 0.0
    
            ∇ = system.∇ .* build_3d_stencil(system, potential.n_datapoints, stencil, system.stencil∇)

        end
    end

    return ∇
end