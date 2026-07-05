function test_3DHarmonicOscillator()

    benchmark_path = base_path * "/benchmark/3DHarmonicOscillator/"

    run_testcase("3DHarmonicOscillator") do
        @suppress Numerov.numerov("input.in")

        compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
        compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
        @test_skip compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat", 1)
    end
end
