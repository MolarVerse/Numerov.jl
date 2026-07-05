function test_1DH2()

    benchmark_path = base_path * "/benchmark/1DH2/"

    run_testcase("1DH2") do
        test_1D_nonreciprocal(benchmark_path)
    end
end
