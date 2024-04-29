function test_1D_nonreciprocal(path::String, benchmark_path::String)
    cd(path)

    cleanup_directory("input.in", "potential.dat")

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
    compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
    compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat", 1)

    cleanup_directory("input.in", "potential.dat")
end