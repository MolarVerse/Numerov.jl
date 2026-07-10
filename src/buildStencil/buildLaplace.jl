"""
	buildΔ(system::System, potential::Potential)

Builds the Laplace operator matrix `Δ` for the given system and potential.

# Arguments

	- `system::System`: The system to build the Laplace operator for.
	- `potential::Potential`: The potential to build the Laplace operator for.

# Returns

	  - `nothing`
"""
function buildΔ(system::System, potential::Potential)

	potential.dimension == 1 && buildΔ_1D(system)
	potential.dimension == 2 && buildΔ_2D(system)
	potential.dimension == 3 && buildΔ_3D(system)

	return nothing
end
