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

    for solver in (Numerov.ARPACK, Numerov.KRYLOV, Numerov.LU)
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

    # the GPU enum value is rejected with a catchable error
    system, output, files = make_solver_structs(Numerov.GPU, n)
    @test_throws ArgumentError Numerov.solveWrapper(system, output, files, H)
end
