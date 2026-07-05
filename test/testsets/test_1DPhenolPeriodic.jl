function test_1DPhenolPeriodic()

    benchmark_path = base_path * "/benchmark/1DPhenolPeriodic/"

    run_testcase("1DPhenolPeriodic") do
        test_1D_nonreciprocal(benchmark_path)
    end
end
