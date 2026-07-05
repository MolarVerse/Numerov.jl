#!/usr/bin/env julia
#
# Command-line launcher for the Numerov solver:
#
#     julia --project=/path/to/Numerov.jl bin/numerov.jl input.in
#
# All output files are written to the current working directory.

using Numerov

if length(ARGS) != 1
    println(stderr, "usage: julia --project=<Numerov.jl dir> bin/numerov.jl <input file>")
    exit(1)
end

numerov(ARGS[1])
