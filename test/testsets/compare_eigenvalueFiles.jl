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