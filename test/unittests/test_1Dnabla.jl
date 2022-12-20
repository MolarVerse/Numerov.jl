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