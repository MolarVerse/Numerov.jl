function unittests()
    @testset "test 1D laplace operator" test_1DΔ()
    @testset "test 1D nabla   operator" test_1D∇()
    @testset "test 2D laplace operator" test_2DΔ()
end

function test_1DΔ()
    system = Numerov.System1D()

    system.stencilΔ = 3
    system.n_datapoints = [10]
    system.periodic = [false]

    Δ = spdiagm( -1 => ones(9),
                  0 => ones(10)*(-2),
                  1 => ones(9))
    
    @test Numerov.buildΔ(system) == Δ

    system.periodic = [true]

    Δ = spdiagm( -9 => ones(1),
                 -1 => ones(9),
                  0 => ones(10)*(-2),
                  1 => ones(9),
                  9 => ones(1),)

    @test Numerov.buildΔ(system) == Δ

end

function test_1D∇()
    system = Numerov.System1D()

    system.stencil∇ = 3
    system.n_datapoints = [10]
    system.periodic = [false]

    ∇ = spdiagm( -1 =>  ones( 9)*(-0.5),
                  0 => zeros(10),
                  1 =>  ones( 9)*( 0.5))
    
    @test Numerov.build∇(system) == ∇

    system.periodic = [true]

    ∇ = spdiagm( -9 =>  ones(1)*( 0.5),
                 -1 =>  ones(9)*(-0.5),
                  0 => zeros(10),
                  1 =>  ones(9)*( 0.5),
                  9 =>  ones(1)*(-0.5),)

    @test Numerov.build∇(system) == ∇

end

function test_2DΔ()
    system = Numerov.System2D()

    system.stencilΔ = 5
    system.n_datapoints = [6,7]
    system.periodic = [false, false]

    Δ = Numerov.buildΔ(system)

    stencil = zeros(system.stencilΔ, system.stencilΔ)
    stencil[:,1] = [ 0.0,  0.0,  -1.0,  0.0,  0.0]
    stencil[:,2] = [ 0.0,  0.0,  16.0,  0.0,  0.0]
    stencil[:,3] = [-1.0, 16.0, -60.0, 16.0, -1.0]
    stencil /= 6.0

    n1 = system.n_datapoints[1]
    n2 = system.n_datapoints[2]

    for j in 0:2
        for i in 0:n1-1-j
            submatrix = Numerov.build_1d_stencil(system, n2, stencil[:,3-j], system.stencilΔ)

            @test Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix
            @test Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix

            Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
            Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
        end
    end

    @test Δ == zeros(n1*n2, n1*n2)

    system.periodic = [true, true]

    Δ = Numerov.buildΔ(system)

    for j in 0:2
        for i in 0:n1-1-j
            submatrix = Numerov.build_1d_stencil(system, n2, stencil[:,3-j], system.stencilΔ)

            @test Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix
            @test Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] == submatrix

            Δ[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
            Δ'[1+i*n2:(i+1)*n2, 1+(i+j)*n2:(i+1+j)*n2] = zeros(n2, n2)
        end
    end

    @test Δ[36:42,1:7]  == Numerov.build_1d_stencil(system, n2, stencil[:,2], system.stencilΔ)
    @test Δ'[36:42,1:7] == Numerov.build_1d_stencil(system, n2, stencil[:,2], system.stencilΔ)
    Δ[36:42,1:7]  = zeros(n2,n2)
    Δ'[36:42,1:7] = zeros(n2, n2)
    
    @test Δ[29:35,1:7]  == Numerov.build_1d_stencil(system, n2, stencil[:,1], system.stencilΔ)
    @test Δ'[29:35,1:7] == Numerov.build_1d_stencil(system, n2, stencil[:,1], system.stencilΔ)
    Δ[29:35,1:7]  = zeros(n2,n2)
    Δ'[29:35,1:7] = zeros(n2,n2)

    @test Δ[36:42,8:14]  == Numerov.build_1d_stencil(system, n2, stencil[:,1], system.stencilΔ)
    @test Δ'[36:42,8:14] == Numerov.build_1d_stencil(system, n2, stencil[:,1], system.stencilΔ)
    Δ[36:42,8:14]  = zeros(n2,n2)
    Δ'[36:42,8:14] = zeros(n2,n2)

    @test Δ == zeros(n1*n2, n1*n2)

end

function test_2D∇()
    system = Numerov.System1D()

    system.stencil∇ = 3
    system.n_datapoints = [10]
    system.periodic = [false]

    ∇ = spdiagm( -1 =>  ones( 9)*(-0.5),
                  0 => zeros(10),
                  1 =>  ones( 9)*( 0.5))
    
    @test Numerov.build∇(system) == ∇

    system.periodic = [true]

    ∇ = spdiagm( -9 =>  ones(1)*( 0.5),
                 -1 =>  ones(9)*(-0.5),
                  0 => zeros(10),
                  1 =>  ones(9)*( 0.5),
                  9 =>  ones(1)*(-0.5),)

    @test Numerov.build∇(system) == ∇

end