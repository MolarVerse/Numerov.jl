function test_2DKronigPenney()

    benchmark_path = base_path * "/benchmark/2DKronigPenney/"

    run_testcase("2DKronigPenney") do
        test_2D_reciprocal(benchmark_path)
        compare_eigenvalueFiles( "bandstructure.dat", benchmark_path * "bandstructure.dat")
    end

    benchmark_path = base_path * "/benchmark/2DKronigPenney_full/"

    run_testcase("2DKronigPenney_full") do
        test_2D_reciprocal(benchmark_path)
    end
end
