####################################################################
#                                                                  #
# basis-invariant subspace checks (GitHub issue #3)                #
#                                                                  #
# Element-wise eigenvector comparisons are ill-defined: degenerate #
# eigenspaces let the solver return any orthonormal basis of the   #
# subspace, and even non-degenerate states carry an arbitrary      #
# sign/phase. These tests compare the spanned subspaces instead,   #
# which is invariant under basis rotation and phase.               #
#                                                                  #
####################################################################

"""
    cluster_by_energy(energies; atol) -> Vector{UnitRange{Int}}

Group the ascending `energies` into clusters of (near-)degenerate levels: a
new cluster starts whenever the gap to the previous energy exceeds `atol`.
Returns the index ranges of the clusters.
"""
function cluster_by_energy(energies::AbstractVector{<:Real}; atol::Real)
    issorted(energies) || throw(ArgumentError("energies must be ascending"))

    clusters = UnitRange{Int}[]
    isempty(energies) && return clusters

    start = 1
    for i in 2:length(energies)
        if energies[i] - energies[i-1] > atol
            push!(clusters, start:(i-1))
            start = i
        end
    end
    push!(clusters, start:length(energies))

    return clusters
end

"""
    orthonormal_basis(A) -> Matrix

Thin QR orthonormalization of the columns of `A`; the returned columns span
the same subspace but are orthonormal in the plain l2 inner product, which
removes the volume-element normalization of the API states.
"""
orthonormal_basis(A::AbstractMatrix) = Matrix(qr(A).Q)[:, 1:size(A, 2)]

"""
    subspace_overlap_deviation(A, B) -> Float64

`opnorm(S'S - I)` with `S = Q1' Q2` built from orthonormalized bases of the
column spans of `A` and `B`. Zero iff the two subspaces coincide; the value
is `max(sin²θᵢ)` over the principal angles θᵢ between the subspaces.
"""
function subspace_overlap_deviation(A::AbstractMatrix, B::AbstractMatrix)
    S = orthonormal_basis(A)' * orthonormal_basis(B)
    return opnorm(S' * S - I)
end

"""
    projector_distance(A, B) -> Float64

Frobenius norm `norm(P1 - P2)` of the difference of the orthogonal
projectors `P = QQ'` onto the column spans of `A` and `B`, evaluated through
the identity `norm(P1 - P2)² = 2m - 2 norm(Q1'Q2)²` so the n×n projectors
are never materialized.
"""
function projector_distance(A::AbstractMatrix, B::AbstractMatrix)
    S = orthonormal_basis(A)' * orthonormal_basis(B)
    m = size(S, 1)
    return sqrt(max(0.0, 2m - 2 * norm(S)^2))
end

"""
    check_cluster_subspaces(sol1, sol2, clusters, volume_element)

For every cluster of (near-)degenerate states assert basis-invariantly that
the two solutions span the same eigenspace:

- the overlap matrix of the orthonormalized cluster bases is unitary,
- the projectors onto the cluster subspaces agree in Frobenius norm,
- every state of one solve lies inside the other solve's cluster eigenspace
  (its projection onto that subspace preserves the norm),
- the columns obey the API normalization `sum(abs2) * volume_element == 1`.
"""
function check_cluster_subspaces(sol1, sol2, clusters, volume_element)
    for r in clusters
        V1 = sol1.states[:, r]
        V2 = sol2.states[:, r]

        for j in axes(V1, 2)
            @test sum(abs2, V1[:, j]) * volume_element ≈ 1.0 atol = 1.0e-8
            @test sum(abs2, V2[:, j]) * volume_element ≈ 1.0 atol = 1.0e-8
        end

        @test subspace_overlap_deviation(V1, V2) < 1.0e-6
        @test projector_distance(V1, V2) < 1.0e-6

        # residual-style physics: re-solving must keep every state of the
        # cluster within the cluster's eigenspace
        Q1 = orthonormal_basis(V1)
        for j in axes(V2, 2)
            v = V2[:, j] / norm(V2[:, j])
            @test norm(Q1' * v) ≈ 1.0 atol = 1.0e-6
        end
    end
end

"""
    test_subspace_invariance()

2D isotropic harmonic oscillator: checks the analytic spectrum
E = nₓ + n_y + 1 → (1, 2, 2, 3, 3, 3) and its degeneracy pattern [1, 2, 3]
via energy clustering, then asserts that the :arpack and :lu solvers agree
per degenerate cluster in the basis-invariant subspace sense. Also
cross-validates the projector-distance identity against explicitly
materialized projectors on this small grid.
"""
function test_subspace_invariance()

    # helper sanity checks on synthetic data
    @test cluster_by_energy([1.0, 2.0, 2.0 + 1.0e-9, 3.0]; atol = 1.0e-6) == [1:1, 2:3, 4:4]
    @test cluster_by_energy(Float64[]; atol = 1.0) == UnitRange{Int}[]
    @test cluster_by_energy([5.0]; atol = 1.0) == [1:1]
    @test cluster_by_energy([1.0, 1.5, 2.1]; atol = 0.7) == [1:3]
    @test_throws ArgumentError cluster_by_energy([2.0, 1.0]; atol = 0.1)

    x = range(-5.0, 5.0; length = 41)
    V = [0.5 * (xi^2 + yi^2) for xi in x, yi in x]

    sol_arpack = solve_schrodinger(V, (x, x); n_eigenvalues = 6, solver = :arpack)
    sol_lu     = solve_schrodinger(V, (x, x); n_eigenvalues = 6, solver = :lu)

    # analytic spectrum E = nₓ + n_y + 1
    for sol in (sol_arpack, sol_lu)
        @test sol.energies ≈ [1.0, 2.0, 2.0, 3.0, 3.0, 3.0] atol = 2.0e-3
        @test sol.kpoint === nothing
        @test eltype(sol.states) <: Real
    end

    # degeneracy pattern [1, 2, 3]
    clusters = cluster_by_energy(sol_arpack.energies; atol = 1.0e-2)
    @test length.(clusters) == [1, 2, 3]
    @test cluster_by_energy(sol_lu.energies; atol = 1.0e-2) == clusters

    dx = step(x)
    check_cluster_subspaces(sol_arpack, sol_lu, clusters, dx^2)

    # on this small grid, cross-check the projector-distance identity
    # against the explicitly materialized projectors P = QQ'
    r  = clusters[end]
    Q1 = orthonormal_basis(sol_arpack.states[:, r])
    Q2 = orthonormal_basis(sol_lu.states[:, r])
    @test projector_distance(sol_arpack.states[:, r], sol_lu.states[:, r]) ≈
          norm(Q1 * Q1' - Q2 * Q2') atol = 1.0e-8
end

"""
    test_subspace_3D()

3D isotropic harmonic oscillator on the 15³ grid of `test_3Dsmoke`: checks
the analytic energies E₀ = 1.5, E₁ = 2.5 and the degeneracy pattern [1, 3],
then asserts :arpack vs :lu subspace agreement per cluster.
"""
function test_subspace_3D()

    x = range(-5.5, 5.5; length = 15)
    V = [0.5 * (xi^2 + yi^2 + zi^2) for xi in x, yi in x, zi in x]

    sol_arpack = solve_schrodinger(V, (x, x, x); n_eigenvalues = 4, solver = :arpack)
    sol_lu     = solve_schrodinger(V, (x, x, x); n_eigenvalues = 4, solver = :lu)

    # analytic spectrum, tolerances as in test_3Dsmoke for the default
    # 9-point stencil
    for sol in (sol_arpack, sol_lu)
        @test sol.energies[1] ≈ 1.5 atol = 0.0032
        for i in 2:4
            @test sol.energies[i] ≈ 2.5 atol = 0.012
        end
    end

    # degeneracy pattern [1, 3]
    clusters = cluster_by_energy(sol_arpack.energies; atol = 1.0e-2)
    @test length.(clusters) == [1, 3]
    @test cluster_by_energy(sol_lu.energies; atol = 1.0e-2) == clusters

    dx = step(x)
    check_cluster_subspaces(sol_arpack, sol_lu, clusters, dx^3)
end

"""
    test_subspace_periodic()

1D cosine band: solves at k = 0 with the :arpack and :krylov solvers and
asserts phase-invariant per-state agreement (|⟨ψₐ|ψ_k⟩| Δx ≈ 1) for the
lowest three non-degenerate states, and that the first row of
`band_structure` reproduces the k = 0 solve.
"""
function test_subspace_periodic()

    n = 64
    a = 2π                                    # lattice constant
    x = range(0.0; step = a / n, length = n)  # periodic grid, endpoint excluded
    V = 2.0 .* (1.0 .- cos.(2π .* x ./ a))

    sol_arpack = solve_schrodinger(V, x; periodic = true, k = [0.0],
                                   n_eigenvalues = 3, solver = :arpack)
    sol_krylov = solve_schrodinger(V, x; periodic = true, k = [0.0],
                                   n_eigenvalues = 3, solver = :krylov)

    @test sol_arpack.kpoint == [0.0]
    @test eltype(sol_arpack.states) <: Complex
    @test sol_arpack.energies ≈ sol_krylov.energies atol = 1.0e-8

    # the lowest three states at k = 0 are non-degenerate
    @test length.(cluster_by_energy(sol_arpack.energies; atol = 1.0e-3)) == [1, 1, 1]

    # the first k-point of the band structure is Γ and its energies
    # reproduce the direct k = 0 solve
    bs = band_structure(V, x; periodic = true, n_kpoints = 3,
                        n_eigenvalues = 3, solver = :arpack)
    @test bs.kpoints[1] ≈ [0.0] atol = 1.0e-12
    @test bs.kpath[1] == 0.0
    @test issorted(bs.kpath)
    @test bs.energies[1, :] ≈ sol_arpack.energies atol = 1.0e-8

    # phase-invariant per-state overlap: |⟨ψₐ|ψ_k⟩| Δx == 1 for states
    # normalized to sum(abs2) * Δx == 1
    dx = step(x)
    for j in eachindex(sol_arpack.energies)
        @test sum(abs2, sol_arpack.states[:, j]) * dx ≈ 1.0 atol = 1.0e-8
        overlap = abs(dot(sol_arpack.states[:, j], sol_krylov.states[:, j])) * dx
        @test overlap ≈ 1.0 atol = 1.0e-6
    end
end
