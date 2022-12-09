function test_1DHarmonicOscillator()

    path           = base_path *  "/testsets/1DHarmonicOscillator/"
    benchmark_path = base_path * "/benchmark/1DHarmonicOscillator/"

    cd(path)

    Numerov.numerov("input.in")

    compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
    compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
    compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat")

    input_files = ["input.in", "potential.dat"]
    rm.(filter(x -> x ∉ input_files, readdir()))

end