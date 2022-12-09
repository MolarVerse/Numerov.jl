include("testsets/1DH2/test_1DH2.jl")

function testsets()
    test_1DH2()
end

function compare_eigenvalueFiles(file1::String, file2::String)

    data1 = readdlm(file1; comments=true)
    data2 = readdlm(file2; comments=true)

    for i in eachindex(data2[1,:])
        tol = 1e-10
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

    @test data1 ≈ data2 atol = 1.0e-10
    
end