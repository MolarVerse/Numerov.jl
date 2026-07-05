"""
	build∇(system::System, potential::Potential)

Builds the gradient operator matrix `∇` for the given system and potential.

# Arguments

	- `system::System`: The system to build the gradient operator for.
	- `potential::Potential`: The potential to build the gradient operator for.

# Returns

	- `nothing`
"""
function build∇(system::System, potential::Potential)

	potential.dimension == 1 && build∇_1D(system)
	potential.dimension == 2 && build∇_2D(system)
	potential.dimension == 3 && build∇_3D(system)

	return nothing
end

"""
	build∇_1D(system::System)

Builds the gradient operator matrix `∇` for the given system in 1D.

# Arguments

	- `system::System`: The system to build the gradient operator for.

# Returns

	- `nothing`
"""
function build∇_1D(system::System) #combine these two functions!

	####################################
	#                                  #
	# retrieve 1d stencil kernel array #
	#                                  #
	####################################

	stencil = get_1d_stencil(system)

	#######################
	#                     #
	# build 1d ∇ operator #
	#                     #
	#######################

	system.∇ =
		build_1d_stencil(system, stencil, system.stencil∇)

	return nothing
end

"""
	build∇_2D(system::System)

Builds the gradient operator matrix `∇` for the given system in 2D.

# Arguments

  - `system::System`: The system to build the gradient operator for.

# Returns

	- `nothing`
"""
function build∇_2D(system::System)

	####################################
	#                                  #
	# retrieve 1d stencil kernel array #
	#                                  #
	####################################

	stencil_1d = get_1d_stencil(system)

	###########################
	#                         #
	# build 2d stencil kernel #
	#                         #
	###########################

	stencil = spzeros(system.stencil∇, system.stencil∇)

	stencil[:, system.stencil∇÷2+1] = stencil_1d
	stencil[system.stencil∇÷2+1, :] = stencil_1d

	#######################
	#                     #
	# build 2d ∇ operator #
	#                     #
	#######################

	system.∇ = build_2d_stencil(
		system,
		system.n_datapoints,
		stencil,
		system.stencil∇,
	)

	return nothing
end

"""
	build∇_3D(system::System)

Builds the gradient operator matrix `∇` for the given system in 3D.

# Arguments

	- `system::System`: The system to build the gradient operator for.

# Returns

	- `nothing`
"""
function build∇_3D(system::System)

	####################################
	#                                  #
	# retrieve 1d stencil kernel array #
	#                                  #
	####################################

	stencil_1d = get_1d_stencil(system)

	###########################
	#                         #
	# build 2d stencil kernel #
	#                         #
	###########################

	stencil =
		zeros(system.stencil∇, system.stencil∇, system.stencil∇)

	stencil[:, system.stencil∇÷2+1, system.stencil∇÷2+1] =
		stencil_1d
	stencil[system.stencil∇÷2+1, :, system.stencil∇÷2+1] =
		stencil_1d
	stencil[system.stencil∇÷2+1, system.stencil∇÷2+1, :] =
		stencil_1d

	#######################
	#                     #
	# build 3d ∇ operator #
	#                     #
	#######################

	system.∇ = build_3d_stencil(
		system,
		system.n_datapoints,
		stencil,
		system.stencil∇,
	)

	return nothing
end


"""
	get_1d_stencil(system::System)

Returns the 1D stencil kernel array for the given system. The stencil kernel array is
determined by the `stencil∇` attribute of the system. Available stencil kernels are
3, 5, 7, 9, and 11 points.

# Arguments

	- `system::System`: The system to get the 1D stencil kernel array for.

# Returns

	- `stencil_1d::Array{Float64, 1}`: The 1D stencil kernel array.
"""
function get_1d_stencil(system::System)

	if system.stencil∇ == 3

		stencil_1d = [-1, 0, 1] ./ 2.0

	elseif system.stencil∇ == 5

		stencil_1d = [1, -8, 0, 8, -1] ./ 12.0

	elseif system.stencil∇ == 7

		stencil_1d = [0, -9, -45, 0, 45, -9, 1] ./ 60.0

	elseif system.stencil∇ == 9

		stencil_1d =
			[3, -32, 168, -672, 0, 672, -168, 32, -3] ./ 840.0

	elseif system.stencil∇ == 11

		stencil_1d =
			[
				-2,
				25,
				-150,
				600,
				-2100,
				0,
				2100,
				-600,
				150,
				-25,
				2,
			] ./ 2520.0

	elseif system.stencil∇ == 13

		throw(ArgumentError("13-point stencil is not supported for the ∇ operator - use 3, 5, 7, 9 or 11!"))

	else

		error("unsupported ∇ stencil size $(system.stencil∇) - use 3, 5, 7, 9 or 11!")

	end

	return stencil_1d
end

"""
	build_1d_stencil(system::System, stencil_1d::Array{Float64, 1}, stencil_size::Int64)

Builds the 1D gradient operator matrix `∇` for the given system.

# Arguments

	- `system::System`: The system to build the gradient operator for.
	- `stencil_1d::Array{Float64, 1}`: The 1D stencil kernel array.
	- `stencil_size::Int64`: The size of the stencil kernel.

# Returns

	- `∇::SparseMatrixCSC{Float64, Int64}`: The 1D gradient operator matrix.
"""
function build∇_k(potential::Potential, system::System, k)
	∇ = spzeros(
		prod(potential.n_datapoints),
		prod(potential.n_datapoints),
	)

	if potential.reciprocal
		if potential.dimension == 1
			∇ = system.∇ * k[1]
		elseif potential.dimension == 2
			if k[2] != 0.0
				∇ = system.∇ * k[2]
				for i in 1:potential.n_datapoints[1]
					∇[
						(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],
						(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],
					] *= k[1] / k[2]
				end
			else
				for i in 1:potential.n_datapoints[1]
					∇[
						(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],
						(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],
					] =
						system.∇[
							(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],
							(i-1)*potential.n_datapoints[2]+1:(i-1)*potential.n_datapoints[2]+potential.n_datapoints[2],
						] * k[1]
				end
			end

		elseif potential.dimension == 3

			stencil = zeros(
				system.stencil,
				system.stencil,
				system.stencil,
			)

			stencil[:, system.stencil÷2+1, system.stencil÷2+1]                  = ones(system.stencil) * k[1]
			stencil[system.stencil÷2+1, :, system.stencil÷2+1]                  = ones(system.stencil) * k[2]
			stencil[system.stencil÷2+1, system.stencil÷2+1, :]                  = ones(system.stencil) * k[3]
			stencil[system.stencil÷2+1, system.stencil÷2+1, system.stencil÷2+1] = 0.0

			∇ =
				system.∇ .* build_3d_stencil(
					system,
					potential.n_datapoints,
					stencil,
					system.stencil∇,
				)

		end
	end

	return ∇
end
