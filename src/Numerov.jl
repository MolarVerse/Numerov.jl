module Numerov

using Unitful
using UnitfulAtomic
using SparseArrays
using LinearAlgebra
using Arpack
using Printf
using PhysicalConstants
using DelimitedFiles

import PhysicalConstants.CODATA2018: h, ħ, N_A, c_0

include("MyUnits.jl")
using .MyUnits

MyUnits.__init__()

abstract type System end

include("datatypes/InputDictionary.jl")
include("datatypes/Potential.jl")
include("datatypes/System.jl")
include("datatypes/Output.jl")

include("readInputFile.jl")
include("checkInput.jl")
include("readPotential.jl")
include("setupSystem.jl")
include("buildStencilMatrices.jl")
include("buildLaplace_1d.jl")
include("buildLaplace_2d.jl")
include("buildLaplace_3d.jl")
include("buildNabla.jl")
include("solve.jl")
include("printResults.jl")
include("main.jl")


end # module Numerov
