function test_2D∇()
    system = Numerov.System()

    system.stencil∇ = 5
    system.n_datapoints = [6,7]
    system.periodic = [false, false]

    potential = Numerov.Potential()
    
    potential.dimension = 2

    Numerov.build∇(system, potential)
    
    ∇ = system.∇

    stencil = zeros(system.stencil∇, system.stencil∇)
    stencil[:,1] = [ 0.0,  0.0,  1,  0.0,  0.0]
    stencil[:,2] = [ 0.0,  0.0,  -8,  0.0,  0.0]
    stencil[:,3] = [1, -8, 0, 8, -1]
    stencil /= 12.0

    n1 = system.n_datapoints[1]
    n2 = system.n_datapoints[2]

    for j in 0:2
        for i in 0:n1-1-j
            submatrix = Numerov.build_1d_stencil(system, stencil[:,3-j], system.stencil∇)
            
            sign1 = 1
            sign2 = -1

            if j != 0
                sign1 = -1
                sign2 = 1
            end

            @test ∇[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == sign1 * submatrix
            @test ∇'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == sign2 * submatrix

            ∇[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
            ∇'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
        end
    end

    @test ∇ == zeros(n1*n2, n1*n2)

    system.periodic = [true, true]

    Numerov.build∇(system, potential)
    
    ∇ = system.∇

    for j in 0:2
        for i in 0:n1-1-j
            submatrix = Numerov.build_1d_stencil(system, stencil[:,3-j], system.stencil∇)

            sign1 = 1
            sign2 = -1

            if j != 0
                sign1 = -1
                sign2 = 1
            end

            @test ∇[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == sign1 * submatrix
            @test ∇'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == sign2 * submatrix

            ∇[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
            ∇'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
        end
    end

    @test ∇[36:42,1:7]  == -Numerov.build_1d_stencil(system, stencil[:,2], system.stencil∇)
    @test ∇'[36:42,1:7] == Numerov.build_1d_stencil(system, stencil[:,2], system.stencil∇)
    ∇[36:42,1:7]  = zeros(n2,n2)
    ∇'[36:42,1:7] = zeros(n2, n2)
    
    @test ∇[29:35,1:7]  == -Numerov.build_1d_stencil(system, stencil[:,1], system.stencil∇)
    @test ∇'[29:35,1:7] == Numerov.build_1d_stencil(system, stencil[:,1], system.stencil∇)
    ∇[29:35,1:7]  = zeros(n2,n2)
    ∇'[29:35,1:7] = zeros(n2,n2)

    @test ∇[36:42,8:14]  == -Numerov.build_1d_stencil(system, stencil[:,1], system.stencil∇)
    @test ∇'[36:42,8:14] == Numerov.build_1d_stencil(system, stencil[:,1], system.stencil∇)
    ∇[36:42,8:14]  = zeros(n2,n2)
    ∇'[36:42,8:14] = zeros(n2,n2)

    @test ∇ == zeros(n1*n2, n1*n2)

end