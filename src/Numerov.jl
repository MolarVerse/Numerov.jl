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

include("datatypes/SolverEnum.jl")
include("datatypes/InputDictionary.jl")
include("datatypes/Potential.jl")
include("datatypes/System.jl")
include("datatypes/Output.jl")
include("datatypes/Files.jl")

include("checkInput/checkFiles.jl")
include("checkInput/checkPotential.jl")
include("checkInput/checkSystem.jl")
include("checkInput/checkOutput.jl")

include("readInputFile.jl")
include("readPotential.jl")
include("setupSystem.jl")
include("setupFiles.jl")
include("buildStencilMatrices.jl")
include("buildLaplace_1d.jl")
include("buildLaplace_2d.jl")
include("buildLaplace_3d.jl")
include("buildNabla.jl")
include("solve.jl")
include("printResults.jl")
include("main.jl")


end # module Numerov
