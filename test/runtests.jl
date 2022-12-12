using Numerov
using Test
using DelimitedFiles
using Suppressor
using Statistics
using SparseArrays

base_path = @__DIR__

include("testsets.jl")
include("unittests.jl")

#testsets()

unittests()