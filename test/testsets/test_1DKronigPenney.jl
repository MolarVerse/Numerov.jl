function test_1DKronigPenney()

    benchmark_path = base_path * "/benchmark/1DKronigPenney/"

    run_testcase("1DKronigPenney") do
        test_1Dbandstructure(benchmark_path)
    end

    benchmark_path = base_path * "/benchmark/1DKronigPenney_23.5u/"

    run_testcase("1DKronigPenney_23.5u") do
        test_1Dbandstructure(benchmark_path)
    end
end
