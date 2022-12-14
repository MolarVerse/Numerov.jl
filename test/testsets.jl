include("testsets/test_1DH2.jl")
include("testsets/test_1DHarmonicOscillator.jl")
include("testsets/test_1DKronigPenney.jl")
include("testsets/test_1DPhenolPeriodic.jl")
include("testsets/test_2DHarmonicOscillator.jl")
include("testsets/test_2DWater.jl")
include("testsets/test_2DKronigPenney.jl")
include("testsets/test_3DHarmonicOscillator.jl")
include("testsets/test_3DKronigPenney.jl")

function testsets()
    # @testset "1D H2" test_1DH2()
    # @testset "1D Harmonic Oscillator" test_1DHarmonicOscillator()
    @testset "1D Kronig Penney" test_1DKronigPenney()
    # @testset "1D Phenol Periodic" test_1DPhenolPeriodic()
    # @testset "2D Harmonic Oscillator" test_2DHarmonicOscillator()
    # @testset "2D Water" test_2DWater()
    @testset "2D Kronig Penney" test_2DKronigPenney()
    @testset "2D Kronig Penney full" test_2DKronigPenney_full() #think of a way to combine this function with the one above
    # @testset "3D Harmonic Oscillator" test_3DHarmonicOscillator()
    # @testset "3D Kronig Penney" test_3DKronigPenney()
end

function compare_eigenvalueFiles(file1::String, file2::String)

    data1 = readdlm(file1; comments=true)
    data2 = readdlm(file2; comments=true)

    for i in eachindex(data2[1,:])
        tol = 1e-9
        if length(split(string(data2[1,i]), ".")[2]) == 6
            tol = 1e-6
        end

        @test data1[:,i] ≈ data2[:,i] atol = tol  
    end

end

function compare_frequenciesFiles(file1::String, file2::String)

    data1 = readdlm(file1; comments=true)
    data2 = readdlm(file2; comments=true)

    data1 = filter(x -> typeof(x) <: Float64, data1)
    data2 = filter(x -> typeof(x) <: Float64, data2)

    @test data1 ≈ data2 rtol = 1.0e-5
    
end

function compare_eigenvectorFiles(file1::String, file2::String, dim::Int64)
    
    data1 = readdlm(file1, comments=true)
    data2 = readdlm(file2, comments=true)

    @test data1[:,1:2] ≈ data2[:,1:2] atol = 1.0e-6

    for i in dim+1:length(data1[1,:])
        @test mean(abs.(data1[:,i])) ≈ mean(abs.(data2[:,i])) atol = 1.0e-6
    end
end