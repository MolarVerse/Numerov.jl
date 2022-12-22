#test for 2d [true,false]-periodic for laplace and nabla and vice versa
function test_2DΔ()
    system = Numerov.System()

    system.stencilΔ = 5
    system.n_datapoints = [6,7]
    system.periodic = [false, false]

    potential = Numerov.Potential()
    
    potential.dimension = 2

    Numerov.buildΔ(system, potential)

    Δ = system.Δ

    stencil = zeros(system.stencilΔ, system.stencilΔ)
    stencil[:,1] = [ 0.0,  0.0,  -1.0,  0.0,  0.0]
    stencil[:,2] = [ 0.0,  0.0,  16.0,  0.0,  0.0]
    stencil[:,3] = [-1.0, 16.0, -60.0, 16.0, -1.0]
    stencil /= 6.0

    n1 = system.n_datapoints[1]
    n2 = system.n_datapoints[2]

    for j in 0:2
        for i in 0:n1-1-j
            submatrix = Numerov.build_1d_stencil(system, stencil[:,3-j], system.stencilΔ)

            @test Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix
            @test Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix

            Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
            Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
        end
    end

    @test Δ == zeros(n1*n2, n1*n2)

    system.periodic = [true, true]

    Numerov.buildΔ(system, potential)

    Δ = system.Δ

    for j in 0:2
        for i in 0:n1-1-j
            submatrix = Numerov.build_1d_stencil(system, stencil[:,3-j], system.stencilΔ)

            @test Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix
            @test Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix

            Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
            Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
        end
    end

    @test Δ[36:42,1:7]  == Numerov.build_1d_stencil(system, stencil[:,2], system.stencilΔ)
    @test Δ'[36:42,1:7] == Numerov.build_1d_stencil(system, stencil[:,2], system.stencilΔ)
    Δ[36:42,1:7]  = zeros(n2,n2)
    Δ'[36:42,1:7] = zeros(n2, n2)
    
    @test Δ[29:35,1:7]  == Numerov.build_1d_stencil(system, stencil[:,1], system.stencilΔ)
    @test Δ'[29:35,1:7] == Numerov.build_1d_stencil(system, stencil[:,1], system.stencilΔ)
    Δ[29:35,1:7]  = zeros(n2,n2)
    Δ'[29:35,1:7] = zeros(n2,n2)

    @test Δ[36:42,8:14]  == Numerov.build_1d_stencil(system, stencil[:,1], system.stencilΔ)
    @test Δ'[36:42,8:14] == Numerov.build_1d_stencil(system, stencil[:,1], system.stencilΔ)
    Δ[36:42,8:14]  = zeros(n2,n2)
    Δ'[36:42,8:14] = zeros(n2,n2)

    @test Δ == zeros(n1*n2, n1*n2)

end