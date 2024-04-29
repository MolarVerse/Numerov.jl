function buildΔ(system::System, potential::Potential)
    potential.dimension == 1 && buildΔ_1D(system)
    potential.dimension == 2 && buildΔ_2D(system)
    potential.dimension == 3 && buildΔ_3D(system)
end