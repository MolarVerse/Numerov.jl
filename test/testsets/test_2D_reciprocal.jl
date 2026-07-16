function test_2D_reciprocal(benchmark_path)
    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")

    files          = filter(x -> startswith(x, "frequencies"), readdir())
    filesBenchmark = filter(x -> startswith(x, "frequencies"), readdir(benchmark_path))

    [compare_frequenciesFiles(files[i], benchmark_path * filesBenchmark[i]) for i in eachindex(files)]
end
