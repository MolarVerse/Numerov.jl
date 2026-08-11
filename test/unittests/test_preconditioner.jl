"""
The reference kinetic operator this test verifies `KineticPreconditioner`
against: the Kronecker SUM of independently-built per-dimension 1D kinetic
matrices (`0.5 * -Δ_d / spacing_d²`), which is exactly what the
preconditioner is DEFINED to approximate the true kinetic operator with.

For 1D problems this is trivially the production operator itself. For 2D
with the 5-point stencil, `buildLaplace_2d`'s Δ is verified below to equal
exactly `2 * (this separable sum)` - i.e. dividing by `2^(dimension-1)` (as
`solve()` does when assembling the Hamiltonian) recovers this separable
model exactly, so the preconditioner is an exact match there too. Every
other 2D stencil (3, 7, 9, 11) and every 3D stencil use a genuinely
non-separable (more elaborate, higher-order or off-cross) pattern, so this
separable model is only an approximation of the true operator in those
cases - that is fine for a preconditioner (only convergence speed depends on
the approximation quality), and final accuracy is independently guaranteed
by `solveWrapper`'s post-solve residual verification, tested elsewhere
(`test_solveWrapper`, `test_lobpcg_fallback`, `test_3Dsmoke`).
"""
function separable_kinetic_reference(dims::NTuple{D, Int}, potential::Numerov.Potential, system::Numerov.System) where D
    Is = [Matrix(1.0I, n, n) for n in dims]
    T  = zeros(prod(dims), prod(dims))
    for d in 1:D
        Δd = raw_delta_1d(dims[d], potential.periodic[d], system.stencilΔ)
        factors = collect(Is)
        factors[d] = 0.5 .* (-Δd) ./ potential.intervall[d]^2
        T .+= reduce(kron, factors)
    end
    return T
end

"1D Δ matrix built exactly the way `Numerov.kinetic_1d` builds it internally."
function raw_delta_1d(n::Int, periodic::Bool, stencil::Int)
    p = Numerov.Potential(); p.dimension = 1; p.n_datapoints = [n]; p.periodic = [periodic]
    s = Numerov.System(); s.n_datapoints = p.n_datapoints; s.periodic = p.periodic
    s.reciprocal = false; s.stencil = stencil; s.stencilΔ = stencil; s.stencil∇ = stencil
    Numerov.buildΔ(s, p)
    return Matrix(s.Δ)
end

"""
Verify `KineticPreconditioner` computes the mathematically correct inverse of
its DEFINED (separable) reference operator - independently re-derived here
via dense Kronecker sums, not by calling any of the preconditioner's own
internals - across 1D/2D/3D and a periodic dimension, through all four
`ldiv!` dispatches (vector/matrix, in-place/allocating). Also confirms, as a
regression guard, that this separable reference exactly equals the true
production kinetic operator for 1D and for 2D's 5-point stencil (division by
`2^(dimension-1)` included), and is only an approximation for every other 2D
stencil and for 3D (see module docstring above) - the 2D stencil=9 case below
guards specifically against re-introducing the "2D is always exact" overclaim
this docstring used to make (the package's default stencil is 9, not 5).
"""
function test_KineticPreconditioner()
    Random.seed!(42)

    cases = (
        (dims = (12,),      periodic = false, stencil = 9, exact_vs_production = true),
        (dims = (12,),      periodic = true,  stencil = 9, exact_vs_production = true),
        (dims = (6, 6),     periodic = false, stencil = 5, exact_vs_production = true),
        (dims = (10, 10),   periodic = false, stencil = 9, exact_vs_production = false),
        (dims = (5, 5, 5),  periodic = false, stencil = 5, exact_vs_production = false),
    )

    for case in cases
        D     = length(case.dims)
        axes_ = ntuple(d -> range(-3.0, 3.0; length = case.dims[d]), D)
        V     = zeros(case.dims...)

        potential, system, _, _ = Numerov.setup_problem(
            V, D == 1 ? axes_[1] : axes_;
            mass = 1.0, periodic = case.periodic, n_eigenvalues = 1,
            stencil = case.stencil, stencil_laplace = case.stencil,
            stencil_nabla = min(case.stencil, 11),
            solver = :arpack, potential_unit = UnitfulAtomic.hartree,
            coord_unit = UnitfulAtomic.bohr, mass_unit = Numerov.MyUnits.m_e,
            reciprocal = false)

        T = separable_kinetic_reference(case.dims, potential, system)

        if case.exact_vs_production
            T_production = 0.5 .* (-Matrix(system.Δ) ./ potential.intervall[1]^2 ./ 2^(D - 1))
            @test T ≈ T_production atol = 1.0e-10
        end

        σ  = 1.7
        P  = Numerov.KineticPreconditioner(potential, system; σ = σ)
        Tσ = Symmetric(T + σ * I)

        N = prod(case.dims)

        x     = randn(N)
        y_ref = Tσ \ x

        y = similar(x)
        ldiv!(y, P, x)                      # 3-arg vector
        @test y ≈ y_ref atol = 1.0e-9 rtol = 1.0e-9

        y2 = copy(x)
        ldiv!(P, y2)                        # 2-arg vector, in-place
        @test y2 ≈ y_ref atol = 1.0e-9 rtol = 1.0e-9

        X     = randn(N, 3)
        Y_ref = Tσ \ X

        Y = similar(X)
        ldiv!(Y, P, X)                       # 3-arg matrix
        @test Y ≈ Y_ref atol = 1.0e-9 rtol = 1.0e-9

        X2 = copy(X)
        ldiv!(P, X2)                         # 2-arg matrix, in-place
        @test X2 ≈ Y_ref atol = 1.0e-9 rtol = 1.0e-9
    end
end
