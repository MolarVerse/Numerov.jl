"""
	build_1d_stencil(system, stencilCoefficients, stencil_size)

Builds a 1D stencil matrix for the given system.

# Arguments

	- `system::System`: The system to build the stencil matrix for.
	- `stencilCoefficients::Array{Float64,1}`: The stencil coefficients.
	- `stencil_size::Int`: The size of the stencil.

# Returns

	- `matrix::SparseMatrixCSC{Float64,Int}`: The 1D stencil matrix.
"""
function build_1d_stencil(
	system,
	stencilCoefficients,
	stencil_size,
)

	n_datapoints = system.n_datapoints[end]

	matrix = spzeros(n_datapoints, n_datapoints)

	for j in 1:stencil_size

		index = j - stencil_size ÷ 2 - 1

		matrix += spdiagm(
			index =>
				ones(n_datapoints - abs(index)) *
				stencilCoefficients[j],
		)

		if system.periodic[end] && index != 0
			matrix += spdiagm(
				-(sign(index) * n_datapoints - index) =>
					ones(abs(index)) * stencilCoefficients[j],
			) #no idea why -ones
		end
	end

	return matrix
end

"""
	accumulate_stencil_blocks!(Is, Js, Vs, block_size, n_outer, total_points,
	                            periodic_outer, sub_matrix_indices, block_triplets)

Shared core of [`build_2d_stencil`](@ref) and [`build_3d_stencil`](@ref):
given the (I, J, V) nonzero triplets of each stencil offset's sub-block, plus
the row/col block bookkeeping the caller works out, append the final
triplets for one "outer axis" level of block placement to `Is`/`Js`/`Vs`.

Building the result via COO triplet accumulation followed by a single
`sparse(...)` call (instead of the equivalent-but-much-slower repeated
`matrix[range_i, range_j] = block` assignment into an already-existing large
sparse matrix - each such assignment can touch and reallocate a large
fraction of the matrix's internal CSC storage) is what makes this fast: it
turns an operation that scaled worse than quadratically in the grid size
into one that is close to linear.

For each row-block `j`, a small per-`j` `Dict` (at most `stencil_size`
entries, discarded every iteration) resolves which stencil offset wins the
placement for a given target column-block *before* any triplets are
appended, exactly reproducing the original's `matrix[...] = block`
last-write-wins overwrite semantics - this matters only when a periodic
wraparound collides with a non-wrapped placement, which can only happen for
grids smaller than the stencil width, but is reproduced correctly regardless.
"""
function accumulate_stencil_blocks!(
	Is,
	Js,
	Vs,
	block_size,
	n_outer,
	total_points,
	periodic_outer,
	sub_matrix_indices,
	block_triplets,
)

	stencil_size = length(sub_matrix_indices)
	winners = Dict{Int, Int}()

	for j in 1:n_outer

		empty!(winners)
		row_offset = (j - 1) * block_size

		for i in 1:stencil_size

			cell_index = j + sub_matrix_indices[i]
			col_start = (cell_index - 1) * block_size + 1

			if cell_index < 1 || cell_index > n_outer
				if periodic_outer
					col_start -= sign(cell_index - 1) * total_points
				else
					continue
				end
			end

			winners[col_start] = i

		end

		for (col_start, i) in winners

			I1, J1, V1 = block_triplets[i]
			col_offset = col_start - 1

			append!(Is, I1 .+ row_offset)
			append!(Js, J1 .+ col_offset)
			append!(Vs, V1)

		end

	end

	return nothing
end

"""
	build_2d_stencil(system, n_datapoints, stencil, stencil_size)

Builds a 2D stencil matrix for the given system.

# Arguments

	- `system::System`: The system to build the stencil matrix for.
	- `n_datapoints::Tuple{Int,Int}`: The number of datapoints in each dimension.
	- `stencil::Array{Array{Float64,1},2}`: The stencil coefficients.
	- `stencil_size::Int`: The size of the stencil.

# Returns

	- `matrix::SparseMatrixCSC{Float64,Int}`: The 2D stencil matrix.
"""
function build_2d_stencil(
	system,
	n_datapoints,
	stencil,
	stencil_size,
)

	total_points = prod(n_datapoints)

	sub_matrix_indices = [i - stencil_size ÷ 2 - 1 for i in 1:stencil_size]
	block_triplets = [
		findnz(build_1d_stencil(system, stencil[i, :], stencil_size))
		for i in 1:stencil_size
	]

	Is = Int[]
	Js = Int[]
	Vs = Float64[]

	accumulate_stencil_blocks!(
		Is,
		Js,
		Vs,
		n_datapoints[2],
		n_datapoints[1],
		total_points,
		system.periodic[end-1],
		sub_matrix_indices,
		block_triplets,
	)

	return sparse(Is, Js, Vs, total_points, total_points)
end

"""
	build_3d_stencil(system, n_datapoints, stencil, stencil_size)

Builds a 3D stencil matrix for the given system.

# Arguments

	- `system::System`: The system to build the stencil matrix for.
	- `n_datapoints::Tuple{Int,Int,Int}`: The number of datapoints in each dimension.
	- `stencil::Array{Array{Array{Float64,1},2},3}`: The stencil coefficients.
	- `stencil_size::Int`: The size of the stencil.

# Returns

	- `matrix::SparseMatrixCSC{Float64,Int}`: The 3D stencil matrix.
"""
function build_3d_stencil(
	system,
	n_datapoints,
	stencil,
	stencil_size,
)

	total_points = prod(n_datapoints)

	sub_matrix_indices = [i - stencil_size ÷ 2 - 1 for i in 1:stencil_size]
	block_triplets = [
		findnz(build_2d_stencil(system, n_datapoints[2:3], stencil[:, :, i], stencil_size))
		for i in 1:stencil_size
	]

	Is = Int[]
	Js = Int[]
	Vs = Float64[]

	accumulate_stencil_blocks!(
		Is,
		Js,
		Vs,
		prod(n_datapoints[2:3]),
		n_datapoints[1],
		total_points,
		system.periodic[end-2],
		sub_matrix_indices,
		block_triplets,
	)

	return sparse(Is, Js, Vs, total_points, total_points)
end
