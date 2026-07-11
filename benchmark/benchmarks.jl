# Benchmark suite for Numerov.jl, consumed by AirspeedVelocity.jl in CI
# (.github/workflows/benchmark.yml) and runnable locally with:
#
#     julia --project=benchmark -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
#     julia --project=benchmark -e 'using BenchmarkTools; include("benchmark/benchmarks.jl"); run(SUITE; verbose=true)'

using BenchmarkTools
using Numerov

const CASES = joinpath(@__DIR__, "..", "test", "testsets")

"""
Run `numerov` on the input files of `test/testsets/<case>/` inside a fresh
temporary directory, with stdout silenced.
"""
function run_case(case::String)
    case_path = joinpath(CASES, case)
    mktempdir() do tmp
        for file in readdir(case_path)
            cp(joinpath(case_path, file), joinpath(tmp, file))
        end
        cd(tmp) do
            redirect_stdout(devnull) do
                numerov("input.in")
            end
        end
    end
end

"""
Write a 3D isotropic harmonic-oscillator case (V = 0.5 r², `n` points per
dimension on [-5.5, 5.5] bohr) into `dir`.
"""
function write_3d_harmonic(dir::String, n::Int)
    xs = range(-5.5, 5.5; length = n)
    open(joinpath(dir, "potential.dat"), "w") do io
        for x in xs, y in xs, z in xs
            println(io, x, " ", y, " ", z, " ", 0.5 * (x^2 + y^2 + z^2))
        end
    end
    open(joinpath(dir, "input.in"), "w") do io
        print(io, """
            potential-file = potential.dat
            potential-unit = hartree
            coord-unit     = bohr
            mass-unit      = me
            reduced-mass   = 1.0
            datapoints     = $n, $n, $n
            n-eigenvalues  = 4
            stencil        = 9
            """)
    end
end

function run_3d_harmonic(n::Int)
    mktempdir() do tmp
        write_3d_harmonic(tmp, n)
        cd(tmp) do
            redirect_stdout(devnull) do
                numerov("input.in")
            end
        end
    end
end

"""
Set up Potential/System structs for the 3D harmonic case without solving,
so operator assembly can be benchmarked in isolation.
"""
function setup_3d_structs(n::Int)
    mktempdir() do tmp
        write_3d_harmonic(tmp, n)
        cd(tmp) do
            potential = Numerov.Potential()
            system    = Numerov.System()
            output    = Numerov.Output()
            files     = Numerov.Files()
            [Numerov.inputDictionary[key] = "" for key in keys(Numerov.inputDictionary)]
            Numerov.readInputFile("input.in")
            Numerov.checkInput(potential)
            Numerov.checkInput(system)
            Numerov.checkInput(files)
            Numerov.checkInput(output)
            Numerov.readPotential(potential, files)
            Numerov.setupSystem(potential, system)
            return potential, system
        end
    end
end

const SUITE = BenchmarkGroup()

SUITE["solve"] = BenchmarkGroup()
SUITE["solve"]["1D_harmonic_201"] = @benchmarkable run_case("1DHarmonicOscillator") seconds = 30 samples = 5
SUITE["solve"]["2D_water"]        = @benchmarkable run_case("2DWater") seconds = 60 samples = 3
SUITE["solve"]["3D_harmonic_15"]  = @benchmarkable run_3d_harmonic(15) seconds = 120 samples = 3

SUITE["bandstructure"] = BenchmarkGroup()
SUITE["bandstructure"]["1D_kronigpenney_10k"] =
    @benchmarkable run_case("1DKronigPenney") seconds = 60 samples = 3

SUITE["assembly"] = BenchmarkGroup()
SUITE["assembly"]["laplacian_3D_15"] = @benchmarkable Numerov.buildΔ(s, p) setup =
    ((p, s) = setup_3d_structs(15)) seconds = 30 samples = 5
