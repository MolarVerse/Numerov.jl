function test_1DΔ()
    system = Numerov.System()

    system.stencilΔ = 3
    system.n_datapoints = [10]
    system.periodic = [false]

    Δ = spdiagm( -1 => ones(9),
                  0 => ones(10)*(-2),
                  1 => ones(9))

    potential = Numerov.Potential()
    
    potential.dimension = 1
    
    Numerov.buildΔ(system, potential)

    @test system.Δ == Δ

    system.periodic = [true]

    Δ = spdiagm( -9 => ones(1),
                 -1 => ones(9),
                  0 => ones(10)*(-2),
                  1 => ones(9),
                  9 => ones(1),)

    Numerov.buildΔ(system, potential)

    @test system.Δ == Δ

end