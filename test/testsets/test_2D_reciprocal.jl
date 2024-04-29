function test_2D_reciprocal(path, benchmark_path)
    cd(path)

    cleanup_directory("input.in", "potential.dat")

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")

    files          = filter(x -> startswith(x, "frequencies"), readdir())
    filesBenchmark = filter(x -> startswith(x, "frequencies"), readdir(benchmark_path))

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")

    files          = filter(x -> startswith(x, "frequencies"), readdir())
    filesBenchmark = filter(x -> startswith(x, "frequencies"), readdir(benchmark_path))

    [compare_frequenciesFiles(files[i], benchmark_path * filesBenchmark[i]) for i in eachindex(files)]

    files          = filter(x -> startswith(x, "eigenvectors_k"), readdir())
    filesBenchmark = filter(x -> startswith(x, "eigenvectors_k"), readdir(benchmark_path))

    [@test_skip compare_eigenvectorFiles(files[i], benchmark_path * filesBenchmark[i], 1) for i in eachindex(files)]

    files          = filter(x -> startswith(x, "imag_eigenvectors_k"), readdir())
    filesBenchmark = filter(x -> startswith(x, "imag_eigenvectors_k"), readdir(benchmark_path))

    [@test_skip compare_eigenvectorFiles(files[i], benchmark_path * filesBenchmark[i], 1) for i in eachindex(files)]
end