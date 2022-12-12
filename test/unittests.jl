function unittests()
    system = Numerov.System1D()

    system.stencilΔ = 3
    system.n_datapoints = [10]
    system.periodic = [false]
    
    @test Numerov.buildΔ(system) == spdiagm(0 => ones(10)*(-2), -1 => ones(9), 1 => ones(9))
end