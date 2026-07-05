function test_2DWater()

    benchmark_path = base_path * "/benchmark/2DWater/"

    run_testcase("2DWater") do
        @suppress Numerov.numerov("input.in")

        compare_eigenvalueFiles( "eigenvalues.dat"         , benchmark_path * "eigenvalues.dat")
        compare_frequenciesFiles("frequencies.dat"         , benchmark_path * "frequencies.dat")
        compare_eigenvectorFiles("eigenvectors.dat"        , benchmark_path * "eigenvectors.dat", 1)
    end
end
