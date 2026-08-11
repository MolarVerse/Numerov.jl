using Numerov
using Test
using Aqua
using DelimitedFiles
using Suppressor
using Statistics
using SparseArrays
using LinearAlgebra
using Random
using Unitful
using UnitfulAtomic

base_path = @__DIR__

include("testsets.jl")
include("unittests.jl")

@testset "Numerov.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(Numerov)
    end

    testsets()

    unittests()
end
