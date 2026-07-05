function compare_eigenvectorFiles(file1::String, file2::String, dim::Int64)
    
    data1 = readdlm(file1, comments=true)
    data2 = readdlm(file2, comments=true)

    @test data1[:,1:2] ≈ data2[:,1:2] atol = 1.0e-6

    for i in dim+1:length(data1[1,:])
        @test mean(abs.(data1[:,i])) ≈ mean(abs.(data2[:,i])) atol = 1.0e-6
    end
end