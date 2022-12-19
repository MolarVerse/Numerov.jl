function test_3DKronigPenney()

    path           = base_path *  "/testsets/3DKronigPenney/"
    benchmark_path = base_path * "/benchmark/3DKronigPenney/"

    cd(path)

    input_files = ["input.in", "potential.dat"]
    rm.(filter(x -> x ∉ input_files, readdir()))

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")

    files          = filter(x -> startswith(x, "frequencies"), readdir())
    filesBenchmark = filter(x -> startswith(x, "frequencies"), readdir(benchmark_path))

    for i in eachindex(files)
        compare_frequenciesFiles(files[i], benchmark_path * filesBenchmark[i])
    end

    files          = filter(x -> startswith(x, "eigenvectors_k"), readdir())
    filesBenchmark = filter(x -> startswith(x, "eigenvectors_k"), readdir(benchmark_path))

    for i in eachindex(files)
        @test_skip compare_eigenvectorFiles(files[i], benchmark_path * filesBenchmark[i], 1)
    end

    files          = filter(x -> startswith(x, "imag_eigenvectors_k"), readdir())
    filesBenchmark = filter(x -> startswith(x, "imag_eigenvectors_k"), readdir(benchmark_path))

    for i in eachindex(files)
        @test_skip compare_eigenvectorFiles(files[i], benchmark_path * filesBenchmark[i], 1)
    end

    compare_eigenvalueFiles( "bandstructure.dat"         , benchmark_path * "bandstructure.dat")

    cleanup_directory("input.in", "potential.dat")
end