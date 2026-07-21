function make_solver_structs(solver, n_eigenvalues::Int)
    system        = Numerov.System()
    system.solver = solver

    output               = Numerov.Output()
    output.n_eigenvalues = n_eigenvalues

    files    = Numerov.Files()
    files.to = Numerov.TimerOutput()

    return system, output, files
end

function check_eigenpairs(H, eigenvalues, eigenvectors, reference, n)
    @test length(eigenvalues) >= n
    @test issorted(real.(eigenvalues))
    @test real.(eigenvalues[1:n]) ≈ reference[1:n] atol = 1.0e-8

    for i in 1:n
        residual = norm(H * eigenvectors[:, i] - eigenvalues[i] * eigenvectors[:, i])
        @test residual < 1.0e-6
    end
end

function test_solveWrapper()
    N = 40
    n = 4

    # real symmetric positive definite: discrete Laplacian + potential diagonal,
    # the Hamiltonian shape produced by non-periodic runs
    Δ = spdiagm(-1 => -ones(N - 1), 0 => 2 * ones(N), 1 => -ones(N - 1))
    V = spdiagm(0 => collect(range(0.1, 2.0; length = N)))
    H = Δ + V

    reference = eigen(Symmetric(Matrix(H))).values

    for solver in (Numerov.ARPACK, Numerov.KRYLOV, Numerov.LOBPCG, Numerov.LU)
        system, output, files = make_solver_structs(solver, n)
        eigenvalues, eigenvectors = Numerov.solveWrapper(system, output, files, H)
        check_eigenpairs(H, eigenvalues, eigenvectors, reference, n)
    end

    # complex Hermitian: the Hamiltonian shape produced by periodic runs
    # (i*S with antisymmetric S is Hermitian)
    S  = spdiagm(1 => 0.05 * ones(N - 1)) - spdiagm(-1 => 0.05 * ones(N - 1))
    Hc = H + im * S
    @test ishermitian(Hc)

    reference_c = eigen(Hermitian(Matrix(Hc))).values

    for solver in (Numerov.ARPACK, Numerov.KRYLOV, Numerov.LU)
        system, output, files = make_solver_structs(solver, n)
        eigenvalues, eigenvectors = Numerov.solveWrapper(system, output, files, Hc)
        check_eigenpairs(Hc, eigenvalues, eigenvectors, reference_c, n)
    end

    # exactly singular matrix (periodic Laplacian, constant null vector):
    # the shift-invert σ sits slightly below zero, so H - σI stays invertible
    # and the zero eigenvalue is still found accurately
    Hp = spdiagm(-1 => -ones(N - 1), 0 => 2 * ones(N), 1 => -ones(N - 1))
    Hp[1, N] = -1.0
    Hp[N, 1] = -1.0

    system, output, files = make_solver_structs(Numerov.ARPACK, n)
    eigenvalues, eigenvectors = Numerov.solveWrapper(system, output, files, Hp)
    reference_p = eigen(Symmetric(Matrix(Hp))).values
    check_eigenpairs(Hp, eigenvalues, eigenvectors, reference_p, n)
    @test abs(real(eigenvalues[1])) < 1.0e-8

    # lobpcg rejects complex Hermitian (periodic) problems with a clear error
    system, output, files = make_solver_structs(Numerov.LOBPCG, n)
    @test_throws ArgumentError Numerov.solveWrapper(system, output, files, Hc)

    # the residual verifier reports machine-precision pairs as tight and
    # corrupted pairs as loose
    vals, vecs = Numerov.solveWrapper(make_solver_structs(Numerov.ARPACK, n)..., H)
    @test Numerov.max_relative_residual(H, vals, vecs, n) < 1.0e-8
    bad = copy(vecs); bad[:, 1] .= randn(N) ./ sqrt(N)
    @test Numerov.max_relative_residual(H, vals, bad, n) > 1.0e-2
    collapsed = copy(vecs); collapsed[:, 1] .= 0.0
    @test Numerov.max_relative_residual(H, vals, collapsed, n) == Inf

    # the GPU enum value is rejected with a catchable error
    system, output, files = make_solver_structs(Numerov.GPU, n)
    @test_throws ArgumentError Numerov.solveWrapper(system, output, files, H)
end

"""
Deterministically exercise both lobpcg safety nets - not by mocking, but by
triggering the real failure modes:

1. IterativeSolvers.lobpcg refuses to run (throws) when the matrix is smaller
   than 3x the requested block size - this reliably fails both retry
   attempts, forcing the arpack fallback.
2. A tiny `lobpcg_maxiter` produces an under-converged (but not thrown)
   result - this is caught by the post-solve residual verifier, which
   re-solves with arpack.

Both must still return results as accurate as calling arpack directly.
"""
function test_lobpcg_fallback()
    Random.seed!(7)

    # (1) too small for the requested block size -> throws on both attempts
    N, n = 12, 1   # nev = n + 5 = 6; N=12 < 3*nev=18 triggers IterativeSolvers'
                   # internal instability guard on every attempt
    Δ = spdiagm(-1 => -ones(N - 1), 0 => 2 * ones(N), 1 => -ones(N - 1))
    V = spdiagm(0 => collect(range(0.1, 2.0; length = N)))
    H = Δ + V
    reference = eigen(Symmetric(Matrix(H))).values

    system, output, files = make_solver_structs(Numerov.LOBPCG, n)
    local eigenvalues, eigenvectors
    @test_logs (:warn, r"lobpcg attempt") match_mode = :any begin
        eigenvalues, eigenvectors = Numerov.solveWrapper(system, output, files, H)
    end
    check_eigenpairs(H, eigenvalues, eigenvectors, reference, n)

    # (2) large enough to run, but capped at 1 lobpcg iteration -> converges
    # nowhere near tol, caught by the residual verifier and re-solved
    N2, n2 = 40, 4
    Δ2 = spdiagm(-1 => -ones(N2 - 1), 0 => 2 * ones(N2), 1 => -ones(N2 - 1))
    V2 = spdiagm(0 => collect(range(0.1, 2.0; length = N2)))
    H2 = Δ2 + V2
    reference2 = eigen(Symmetric(Matrix(H2))).values

    system2, output2, files2 = make_solver_structs(Numerov.LOBPCG, n2)
    local eigenvalues2, eigenvectors2
    @test_logs (:warn, r"exceed the residual tolerance") match_mode = :any begin
        eigenvalues2, eigenvectors2 = Numerov.solveWrapper(system2, output2, files2, H2; lobpcg_maxiter = 1)
    end
    check_eigenpairs(H2, eigenvalues2, eigenvectors2, reference2, n2)
end

"""
Non-lobpcg solvers have no fallback to escalate to (arpack is already the top
of the ladder), so an under-converged Krylov result only warns rather than
re-solving. Deterministically trigger this with a tiny `krylov_maxiter`
rather than hoping a normal run happens to under-converge.
"""
function test_residual_warning()
    N, n = 60, 4
    Δ = spdiagm(-1 => -ones(N - 1), 0 => 2 * ones(N), 1 => -ones(N - 1))
    V = spdiagm(0 => collect(range(0.1, 2.0; length = N)))
    H = Δ + V

    system, output, files = make_solver_structs(Numerov.KRYLOV, n)
    local eigenvalues, eigenvectors
    @test_logs (:warn, r"residuals are larger than expected") match_mode = :any begin
        eigenvalues, eigenvectors = Numerov.solveWrapper(system, output, files, H; krylov_maxiter = 1)
    end
    # this is a deliberately broken scenario purely to exercise the warning
    # path - just confirm the residual it reports is indeed loose
    @test Numerov.max_relative_residual(H, eigenvalues, eigenvectors, n) > 1.0e-6
end

"""
The lobpcg->arpack escalation (triggered by a loose residual, not by an
exception) must itself be re-verified rather than returned on trust. Force
lobpcg to under-converge (`lobpcg_maxiter = 1`, as in `test_lobpcg_fallback`)
on a Hamiltonian whose extreme diagonal dynamic range also makes the arpack
rescue itself land above the residual tolerance - `solve_arpack`'s shift
heuristic scales with the largest diagonal entry, so a single huge spike
degrades shift-invert accuracy for the low-lying eigenpairs actually wanted.
This must produce a second, distinct warning naming the rescue result
itself as suspect, not silence.
"""
function test_lobpcg_arpack_rescue_reverified()
    Random.seed!(3)

    N, n = 40, 4
    Δ = spdiagm(-1 => -ones(N - 1), 0 => 2 * ones(N), 1 => -ones(N - 1))
    V = collect(range(0.1, 2.0; length = N))
    V[20] = 1.0e10
    H = Δ + spdiagm(0 => V)

    system, output, files = make_solver_structs(Numerov.LOBPCG, n)
    local eigenvalues, eigenvectors
    @test_logs (:warn, r"re-solving with arpack") (:warn, r"arpack rescue itself exceeds") match_mode = :any begin
        eigenvalues, eigenvectors = Numerov.solveWrapper(system, output, files, H; lobpcg_maxiter = 1)
    end

    # the rescue result is still returned (nothing better to fall back to) -
    # confirm it is indeed the loose result the warning describes, not a
    # spuriously-triggered warning on an otherwise-fine result
    @test Numerov.max_relative_residual(H, eigenvalues, eigenvectors, n) > 1.0e-6
end
