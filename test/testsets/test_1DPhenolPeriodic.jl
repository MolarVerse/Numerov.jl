function test_1DPhenolPeriodic()

    path           = base_path *  "/testsets/1DPhenolPeriodic/"
    benchmark_path = base_path * "/benchmark/1DPhenolPeriodic/"

    cd(path)

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
    compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
    compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat", 1)

    input_files = ["input.in", "potential.dat"]
    rm.(filter(x -> x ∉ input_files, readdir()))

end