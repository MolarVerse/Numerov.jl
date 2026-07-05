include("testsets/compare_eigenvalueFiles.jl")
include("testsets/compare_eigenvectorFiles.jl")
include("testsets/compare_frequenciesFiles.jl")

include("testsets/test_1D_nonreciprocal.jl")
include("testsets/test_1D_bandstructure.jl")
include("testsets/test_2D_reciprocal.jl")

include("testsets/test_1DH2.jl")
include("testsets/test_1DHarmonicOscillator.jl")
include("testsets/test_1DKronigPenney.jl")
include("testsets/test_1DPhenolPeriodic.jl")
include("testsets/test_2DHarmonicOscillator.jl")
include("testsets/test_2DWater.jl")
include("testsets/test_2DKronigPenney.jl")
include("testsets/test_3DHarmonicOscillator.jl")
include("testsets/test_3DKronigPenney.jl")
include("testsets/test_3Dsmoke.jl")
include("testsets/test_solverVariants.jl")

"""
    run_testcase(f, case)

Copy all input files of `test/testsets/<case>/` into a fresh temporary
directory and run `f` inside it, so the solver output never ends up in
the package directory.
"""
function run_testcase(f::Function, case::String)
    case_path = joinpath(base_path, "testsets", case)
    mktempdir() do tmp
        for file in readdir(case_path)
            cp(joinpath(case_path, file), joinpath(tmp, file))
        end
        cd(f, tmp)
    end
end

function testsets()
    @testset "1D H2" test_1DH2()
    @testset "1D Harmonic Oscillator" test_1DHarmonicOscillator()
    @testset "1D Kronig Penney" test_1DKronigPenney()
    @testset "1D Phenol Periodic" test_1DPhenolPeriodic()
    @testset "2D Harmonic Oscillator" test_2DHarmonicOscillator()
    @testset "2D Water" test_2DWater()
    @testset "2D Kronig Penney" test_2DKronigPenney()
    @testset "3D smoke (analytic harmonic oscillator)" test_3Dsmoke()
    @testset "Solver variants" test_solverVariants()

    if get(ENV, "NUMEROV_TEST_FULL", "") == "true"
        @testset "3D Harmonic Oscillator" test_3DHarmonicOscillator()
        @testset "3D Kronig Penney" test_3DKronigPenney()
    else
        @info "Skipping slow testsets \"3D Harmonic Oscillator\" and \"3D Kronig Penney\"; set ENV[\"NUMEROV_TEST_FULL\"] = \"true\" to run them."
    end
end
