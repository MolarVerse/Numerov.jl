function test_1DHarmonicOscillator()

    benchmark_path = base_path * "/benchmark/1DHarmonicOscillator/"

    run_testcase("1DHarmonicOscillator") do
        test_1D_nonreciprocal(benchmark_path)
    end
end
