"""
    run_solver_variant(f, case_path, solver)

Copy all input files of `case_path` into a fresh temporary directory,
append a `solver = <solver>` entry to the copied `input.in` and run `f`
inside it, so the solver output never ends up in the package directory.
"""
function run_solver_variant(f::Function, case_path::String, solver::String)
    mktempdir() do tmp
        for file in readdir(case_path)
            cp(joinpath(case_path, file), joinpath(tmp, file))
        end
        open(joinpath(tmp, "input.in"), "a") do io
            println(io, "solver = ", solver)
        end
        cd(f, tmp)
    end
end

function test_solverVariants()

    benchmark_path = base_path * "/benchmark/1DHarmonicOscillator/"
    case_path      = joinpath(base_path, "testsets", "1DHarmonicOscillator")

    reference = vec(readdlm(benchmark_path * "eigenvalues.dat"; comments=true))

    # non-default eigensolver branches of solveWrapper (src/solve.jl);
    # the default branch (empty solver entry -> arpack) is covered by the
    # regular testsets
    for solver in ["krylov", "lu"]
        run_solver_variant(case_path, solver) do
            @suppress Numerov.numerov("input.in")

            @test isfile("eigenvalues.dat")
            eigenvalues = vec(readdlm("eigenvalues.dat"; comments=true))

            @test length(eigenvalues) == length(reference)
            # krylov and lu reproduce the arpack benchmark to ~1e-13 for this
            # system (measured); 1e-6 leaves margin for other platforms
            @test eigenvalues ≈ reference atol = 1e-6
        end
    end

    # the cuda solver is rejected during input validation (checkSolver in
    # src/checkInput/checkSystem.jl), so solveWrapper is never reached
    run_solver_variant(case_path, "cuda") do
        @test_throws ArgumentError @suppress Numerov.numerov("input.in")
    end
end
