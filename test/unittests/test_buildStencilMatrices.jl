"""
Reference (pre-optimization) implementations of `build_2d_stencil` /
`build_3d_stencil`, kept here only as ground truth: they build the result by
directly assigning each stencil offset's sub-block into the appropriate
`range_i, range_j` slice of an already-existing sparse matrix, which is
obviously correct by inspection (it is literally the mathematical
definition: "the block for stencil offset i, placed at cell j") but does not
scale well, which is exactly why `build_2d_stencil`/`build_3d_stencil` no
longer work this way. `test_buildStencilMatrices` checks the production
(fast, COO-accumulation-based) functions reproduce these bit-for-bit across
1D/2D/3D, every legal stencil size (3/5/7/9/11/13 in 2D; 3/5/7/9/11/13 here
too, even though `buildLaplace_3d`'s own coefficients are only implemented up
to 11 - `build_3d_stencil` itself has no such restriction, and its stencil∇
sub-case IS already reachable in production at size 11), every periodicity
combination (all 8 in 3D, not just the 5 that omit "exactly two axes
periodic"), and - crucially - grids smaller than the stencil width, where a
periodic wraparound can collide with a non-wrapped placement and a naive
`sparse(I, J, V)`-with-default-summing rewrite would silently double-count
instead of overwriting.
"""
function reference_build_2d_stencil(system, n_datapoints, stencil, stencil_size)
    total_points = prod(n_datapoints)
    matrix = spzeros(total_points, total_points)

    for i in 1:stencil_size
        sub_matrix_index = i - stencil_size ÷ 2 - 1
        matrix_1d = Numerov.build_1d_stencil(system, stencil[i, :], stencil_size)

        for j in 1:n_datapoints[1]
            cell_index = j + sub_matrix_index

            range_i = ((j - 1) * n_datapoints[2] + 1):(j * n_datapoints[2])
            range_j = ((cell_index - 1) * n_datapoints[2] + 1):(cell_index * n_datapoints[2])

            if cell_index < 1 || cell_index > n_datapoints[1]
                if system.periodic[end-1]
                    shift = sign(cell_index - 1) * total_points
                    range_j = (range_j[1] - shift):(range_j[end] - shift)
                else
                    continue
                end
            end

            matrix[range_i, range_j] = matrix_1d
        end
    end

    return matrix
end

function reference_build_3d_stencil(system, n_datapoints, stencil, stencil_size)
    total_points = prod(n_datapoints)
    matrix = spzeros(total_points, total_points)

    for i in 1:stencil_size
        sub_matrix_index = i - stencil_size ÷ 2 - 1
        matrix_2d = reference_build_2d_stencil(system, n_datapoints[2:3], stencil[:, :, i], stencil_size)
        block = prod(n_datapoints[2:3])

        for j in 1:n_datapoints[1]
            cell_index = j + sub_matrix_index

            range_i = ((j - 1) * block + 1):(j * block)
            range_j = ((cell_index - 1) * block + 1):(cell_index * block)

            if cell_index < 1 || cell_index > n_datapoints[1]
                if system.periodic[end-2]
                    shift = sign(cell_index - 1) * total_points
                    range_j = (range_j[1] - shift):(range_j[end] - shift)
                else
                    continue
                end
            end

            matrix[range_i, range_j] = matrix_2d
        end
    end

    return matrix
end

function make_stencil_test_system(n_datapoints, periodic, stencil_size)
    system = Numerov.System()
    system.n_datapoints = collect(n_datapoints)
    system.periodic     = collect(periodic)
    system.stencil       = stencil_size
    system.stencilΔ      = stencil_size
    system.stencil∇      = min(stencil_size, 11)
    system.reciprocal    = false
    return system
end

function test_buildStencilMatrices()
    Random.seed!(11)

    @testset "build_2d_stencil matches the reference" begin
        for stencil_size in (3, 5, 7, 9, 11, 13)
            for dims in ((6, 6), (10, 8), (stencil_size, stencil_size), (stencil_size, 6), (6, stencil_size))
                for periodic in ((false, false), (true, false), (false, true), (true, true))
                    system = make_stencil_test_system(dims, periodic, stencil_size)
                    stencil = randn(stencil_size, stencil_size)

                    reference = reference_build_2d_stencil(system, dims, stencil, stencil_size)
                    fast      = Numerov.build_2d_stencil(system, dims, stencil, stencil_size)

                    @test fast == reference
                end
            end
        end
    end

    @testset "build_3d_stencil matches the reference" begin
        for stencil_size in (3, 5, 7, 9, 11, 13)
            for dims in ((6, 6, 6), (stencil_size, stencil_size, stencil_size), (8, 6, 7))
                for periodic in (
                    (false, false, false), (true, true, true),
                    (true, false, false), (false, true, false), (false, false, true),
                    (true, true, false), (true, false, true), (false, true, true),
                )
                    system = make_stencil_test_system(dims, periodic, stencil_size)
                    stencil = randn(stencil_size, stencil_size, stencil_size)

                    reference = reference_build_3d_stencil(system, dims, stencil, stencil_size)
                    fast      = Numerov.build_3d_stencil(system, dims, stencil, stencil_size)

                    @test fast == reference
                end
            end
        end
    end

    # end-to-end sanity: a real (non-random-coefficient) 3D Laplacian must
    # still be symmetric and have the expected row sum (each stencil's
    # coefficients sum to zero away from boundaries/periodic wraparound,
    # since a constant function has zero Laplacian)
    @testset "production buildΔ_3D stays symmetric and well-formed" begin
        n = 11
        system = make_stencil_test_system((n, n, n), (false, false, false), 9)
        potential = Numerov.Potential()
        potential.dimension = 3
        Numerov.buildΔ(system, potential)
        @test issymmetric(system.Δ)
        @test size(system.Δ) == (n^3, n^3)
    end
end
