function test_2DHarmonicOscillator()

    path           = base_path *  "/testsets/2DHarmonicOscillator/"
    benchmark_path = base_path * "/benchmark/2DHarmonicOscillator/"

    cd(path)

    @suppress Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
    compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
    @test_skip compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat", 1)

    input_files = ["input.in", "potential.dat"]
    rm.(filter(x -> x ∉ input_files, readdir()))

end