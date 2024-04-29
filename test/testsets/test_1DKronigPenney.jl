function test_1DKronigPenney()

    path           = base_path *  "/testsets/1DKronigPenney/"
    benchmark_path = base_path * "/benchmark/1DKronigPenney/"

    test_1Dbandstructure(path, benchmark_path)

    path           = base_path *  "/testsets/1DKronigPenney_23.5u/"
    benchmark_path = base_path * "/benchmark/1DKronigPenney_23.5u/"

    test_1Dbandstructure(path, benchmark_path)
end