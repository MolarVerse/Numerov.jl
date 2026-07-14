function test_3DKronigPenney()

    benchmark_path = base_path * "/benchmark/3DKronigPenney/"

    run_testcase("3DKronigPenney") do
        @suppress Numerov.numerov("input.in")

        compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")

        files          = filter(x -> startswith(x, "frequencies"), readdir())
        filesBenchmark = filter(x -> startswith(x, "frequencies"), readdir(benchmark_path))

        [compare_frequenciesFiles(files[i], benchmark_path * filesBenchmark[i]) for i in eachindex(files)]

        compare_eigenvalueFiles( "bandstructure.dat"         , benchmark_path * "bandstructure.dat")
    end
end
