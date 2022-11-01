module Numerov

using Unitful
using UnitfulAtomic
using SparseArrays
using LinearAlgebra
using Arpack
using Printf
using PhysicalConstants

import PhysicalConstants.CODATA2018: h, ħ, N_A, c_0

include("MyUnits.jl")
using .MyUnits

abstract type System end

include("datatypes/InputDictionary.jl")
include("datatypes/Potential.jl")
include("datatypes/System.jl")
include("datatypes/Output.jl")

include("readInputFile.jl")
include("checkInput.jl")
include("readPotential.jl")
include("setupSystem.jl")
include("buildLaplace.jl")
include("solve.jl")
include("printResults.jl")
include("main.jl")


end # module Numerov
