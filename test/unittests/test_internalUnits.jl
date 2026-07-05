function test_internalUnits()
    potential = Numerov.Potential()
    
    Numerov.set_internalUnits(potential)
 
    @test potential.internalElemEnergy == u"hartree"
    @test potential.internalElemCoords == u"bohr"
    @test potential.internalElemMass   == u"m_e"
 end