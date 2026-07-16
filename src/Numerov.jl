module Numerov
    
    using Comonicon
    using Unitful
    using UnitfulAtomic
    using SparseArrays
    using LinearAlgebra
    using Arpack
    using Printf
    using PhysicalConstants
    using DelimitedFiles
    using TimerOutputs
    using KrylovKit
    using StatsBase

    import PhysicalConstants.CODATA2018: h, ħ, N_A, c_0

    export numerov
    export solve_schrodinger, band_structure
    export SchrodingerSolution, BandStructure

    include("MyUnits.jl")
    using .MyUnits

    MyUnits.__init__()

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

    include("buildStencil/buildStencilMatrices.jl")
    include("buildStencil/buildLaplace_1d.jl")
    include("buildStencil/buildLaplace_2d.jl")
    include("buildStencil/buildLaplace_3d.jl")
    include("buildStencil/buildLaplace.jl")
    include("buildStencil/buildNabla.jl")

    include("printResults/printEigenvalues.jl")
    include("printResults/printEigenvectors.jl")
    include("printResults/printFrequencies.jl")
    include("printResults/printBandStructure.jl")

    include("printLogFile/init_logfile.jl")
    include("printLogFile/sparseInfo.jl")
    include("printLogFile/fileInfo.jl")
    include("printLogFile/systemInfo.jl")

    include("k_paths.jl")
    include("readInputFile.jl")
    include("readPotential.jl")

    include("setupSystem.jl")

    include("normalize.jl")
    include("solve.jl")

    include("main.jl")
    include("api.jl")


end # module Numerov
