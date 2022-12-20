function test_2DKronigPenney()

    path           = base_path *  "/testsets/2DKronigPenney/"
    benchmark_path = base_path * "/benchmark/2DKronigPenney/"

    test_2D_reciprocal(path, benchmark_path)
    compare_eigenvalueFiles( "bandstructure.dat", benchmark_path * "bandstructure.dat")
    cleanup_directory("input.in", "potential.dat")

    path           = base_path *  "/testsets/2DKronigPenney_full/"
    benchmark_path = base_path * "/benchmark/2DKronigPenney_full/"

    test_2D_reciprocal(path, benchmark_path)
    cleanup_directory("input.in", "potential.dat")
end