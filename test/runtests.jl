using Numerov
using Test
using DelimitedFiles

base_path = @__DIR__

include("testsets.jl")
include("test_1DH2.jl")
include("test_1DHarmonicOscillator.jl")

testsets()