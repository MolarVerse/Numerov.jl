function test_3DHarmonicOscillator()

    path           = base_path *  "/testsets/3DHarmonicOscillator/"
    benchmark_path = base_path * "/benchmark/3DHarmonicOscillator/"

    cd(path)

    input_files = ["input.in", "potential.dat"]
    rm.(filter(x -> x ∉ input_files, readdir()))

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
    compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
    @test_skip compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat", 1)

    cleanup_directory("input.in", "potential.dat")
end