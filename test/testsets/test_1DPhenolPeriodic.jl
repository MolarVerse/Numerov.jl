function test_1DPhenolPeriodic()

    path           = base_path *  "/testsets/1DPhenolPeriodic/"
    benchmark_path = base_path * "/benchmark/1DPhenolPeriodic/"

    test_1D_nonreciprocal(path, benchmark_path)
end