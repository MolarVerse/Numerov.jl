using Numerov
using Test
using DelimitedFiles
using Suppressor
using Statistics
using SparseArrays
using Unitful
using UnitfulAtomic

base_path = @__DIR__

include("testsets.jl")
include("unittests.jl")

# testsets()

unittests()