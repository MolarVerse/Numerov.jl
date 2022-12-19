function test_1DHarmonicOscillator()

    path           = base_path *  "/testsets/1DHarmonicOscillator/"
    benchmark_path = base_path * "/benchmark/1DHarmonicOscillator/"

    test_1D_nonreciprocal(path, benchmark_path)
end