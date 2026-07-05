function compare_frequenciesFiles(file1::String, file2::String)

    data1 = readdlm(file1; comments=true)
    data2 = readdlm(file2; comments=true)

    data1 = filter(x -> typeof(x) <: Float64, data1)
    data2 = filter(x -> typeof(x) <: Float64, data2)

    @test data1 ≈ data2 rtol = 1.0e-5
end