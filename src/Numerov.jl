module Numerov

using Unitful
using UnitfulAtomic
using SparseArrays
using LinearAlgebra

abstract type System end

include("datatypes/InputDictionary.jl")
include("datatypes/Potential.jl")
include("datatypes/System.jl")

include("readInputFile.jl")
include("checkInput.jl")
include("readPotential.jl")
include("setupSystem.jl")
include("buildLaplace.jl")
include("main.jl")


end # module Numerov
