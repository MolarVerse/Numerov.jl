##############################################################################
#                                                                            #
# unit tests for the input file parsing (src/readInputFile.jl) and the      #
# input validation (src/checkInput/*.jl)                                     #
#                                                                            #
##############################################################################

using Numerov.MyUnits # for u"m_e" and u"u" in the mass unit tests

# reset the global input dictionary to its default (empty) values the same
# way numerov() does at the start of src/main.jl
function reset_inputDictionary()
    for key in keys(Numerov.inputDictionary)
        Numerov.inputDictionary[key] = ""
    end
end

# write `content` to an input file inside a fresh temporary directory and
# parse it with Numerov.readInputFile after resetting the global dictionary
function parse_inputString(content::String)
    reset_inputDictionary()
    mktempdir() do tmp
        inputFileName = joinpath(tmp, "input.in")
        write(inputFileName, content)
        Numerov.readInputFile(inputFileName)
    end
end

# parse `content` and run the check of the potential related keywords on it
function checkPotential_from_input(content::String)
    parse_inputString(content)
    potential = Numerov.Potential()
    Numerov.checkInput(potential)
    return potential
end

# parse `content` and run the check of the system related keywords on it
function checkSystem_from_input(content::String)
    parse_inputString(content)
    system = Numerov.System()
    Numerov.checkInput(system)
    return system
end

# parse `content` and run the check of the output related keywords on it
function checkOutput_from_input(content::String)
    parse_inputString(content)
    output = Numerov.Output()
    Numerov.checkInput(output)
    return output
end

# run a full numerov() calculation on the input file `content` inside a
# fresh temporary directory so no output files end up in the repository
function run_numerov_on_input(content::String; potentialFileContent::Union{Nothing, String}=nothing)
    mktempdir() do tmp
        cd(tmp) do
            write("input.in", content)
            potentialFileContent === nothing || write("potential.dat", potentialFileContent)
            @suppress Numerov.numerov("input.in")
        end
    end
end

# tiny 1D harmonic potential on a 20 point grid in "x V" format
function harmonic_potential_1D()
    return join(["$(x) $((x - 1.0)^2)" for x in 0.0:0.1:1.9], "\n") * "\n"
end

function test_readInputFile_errors()
    # unknown keyword
    @test_throws ArgumentError parse_inputString("some-keyword = 1\n")

    # keyword defined multiple times - also case insensitively
    @test_throws ArgumentError parse_inputString("stencil = 3\nstencil = 5\n")
    @test_throws ArgumentError parse_inputString("Stencil = 3\nstencil = 5\n")

    # second token of a line has to be a "="
    @test_throws ArgumentError parse_inputString("stencil : 3\n")

    # too few tokens in a line
    @test_throws ArgumentError parse_inputString("stencil =\n")
    @test_throws ArgumentError parse_inputString("stencil\n")
end

function test_readInputFile_valid()
    # keywords are case insensitive, comments and blank lines are ignored and
    # multi token values are joined again with single spaces
    parse_inputString("STENCIL = 3 # a comment\n\npotential-file = pot.dat\nreduced-mass = 1.0, 2.0\n")

    @test Numerov.inputDictionary["stencil"]        == "3"
    @test Numerov.inputDictionary["potential-file"] == "pot.dat"
    @test Numerov.inputDictionary["reduced-mass"]   == "1.0, 2.0"
end

function test_checkPotential_errors()
    # unknown unit values
    @test_throws ArgumentError checkPotential_from_input("potential-unit = joule\n")
    @test_throws ArgumentError checkPotential_from_input("coord-unit = meter\n")
    @test_throws ArgumentError checkPotential_from_input("mass-unit = kg\n")

    # unknown band-structure option
    @test_throws ArgumentError checkPotential_from_input("band-structure = maybe\n")

    # number of k-points has to be > 1
    @test_throws ArgumentError checkPotential_from_input("k-points = 1\n")

    # "," parses to an empty vector for reduced-mass and datapoints
    @test_throws ArgumentError checkPotential_from_input("reduced-mass = ,\n")
    @test_throws ArgumentError checkPotential_from_input("datapoints = ,\n")

    # reading k-points from a file is only allowed if no k-points are given
    @test_throws ArgumentError checkPotential_from_input("read-k-points = true\nk-points = 5\nk-points-file = k.dat\n")

    # reading k-points from a file requires the k-points-file keyword
    @test_throws ArgumentError checkPotential_from_input("read-k-points = true\n")
end

function test_checkPotential_accepted()
    # "off" and "false" are accepted and disable the band structure
    for value in ["off", "false", "OFF"]
        potential = checkPotential_from_input("band-structure = $(value)\n")
        @test potential.bandStructure == false
    end

    # "on" and "true" are accepted and enable the band structure
    for value in ["on", "true"]
        potential = checkPotential_from_input("band-structure = $(value)\n")
        @test potential.bandStructure == true
    end

    # "me" and "m_e" both stand for the electron mass unit
    for value in ["me", "m_e"]
        potential = checkPotential_from_input("mass-unit = $(value)\n")
        @test potential.massUnit == u"m_e"
    end

    # "unit" and "g/mol" both map to unified atomic mass units
    for value in ["unit", "g/mol"]
        potential = checkPotential_from_input("mass-unit = $(value)\n")
        @test potential.massUnit == u"u"
    end
end

function test_checkSystem_errors()
    # stencil and stencil-laplace have to be one of 3, 5, 7, 9, 11, 13
    @test_throws ArgumentError checkSystem_from_input("stencil = 4\n")
    @test_throws ArgumentError checkSystem_from_input("stencil-laplace = 4\n")

    # stencil-nabla has to be one of 3, 5, 7, 9, 11
    @test_throws ArgumentError checkSystem_from_input("stencil-nabla = 13\n")

    # the cuda solver is recognised but not implemented
    @test_throws ArgumentError checkSystem_from_input("solver = cuda\n")

    # unknown solver
    @test_throws ArgumentError checkSystem_from_input("solver = magic\n")
end

function test_checkSystem_accepted()
    # defaults
    @test checkSystem_from_input("").stencil == 9
    @test checkSystem_from_input("").solver  == Numerov.ARPACK

    # solver values are case insensitive
    @test checkSystem_from_input("solver = ARPACK\n").solver == Numerov.ARPACK
    @test checkSystem_from_input("solver = krylov\n").solver == Numerov.KRYLOV
    @test checkSystem_from_input("solver = lu\n").solver     == Numerov.LU

    @test checkSystem_from_input("stencil = 3\n").stencil        == 3
    @test checkSystem_from_input("stencil-nabla = 11\n").stencil∇ == 11
end

function test_checkOutput()
    # number of eigenvalues has to be >= 1
    @test_throws ArgumentError checkOutput_from_input("n-eigenvalues = 0\n")
    @test_throws ArgumentError checkOutput_from_input("n-eigenvalues = -3\n")

    # default and explicit values
    @test checkOutput_from_input("").n_eigenvalues                    == 5
    @test checkOutput_from_input("n-eigenvalues = 7\n").n_eigenvalues == 7
end

function test_numerov_inputValidation()
    # parsing errors are also raised by a full numerov() run
    @test_throws ArgumentError run_numerov_on_input("some-keyword = 1\n")

    # missing potential-file keyword (checked after the potential and system input)
    @test_throws ArgumentError run_numerov_on_input("n-eigenvalues = 3\n")

    # invalid values are caught by the checks before the potential file is accessed
    @test_throws ArgumentError run_numerov_on_input("potential-file = potential.dat\npotential-unit = joule\n")
    @test_throws ArgumentError run_numerov_on_input("potential-file = potential.dat\nsolver = cuda\n")
    @test_throws ArgumentError run_numerov_on_input("potential-file = potential.dat\nstencil = 4\n")
    @test_throws ArgumentError run_numerov_on_input("potential-file = potential.dat\nn-eigenvalues = 0\n")

    # a nonexistent potential file is only detected when it is opened
    @test_throws SystemError run_numerov_on_input("potential-file = does-not-exist.dat\n")

    # band-structure = false is accepted - the full calculation runs through
    @test begin
        run_numerov_on_input("""
                             potential-file = potential.dat
                             band-structure = false
                             solver         = lu
                             n-eigenvalues  = 2
                             mass-unit      = m_e
                             """; potentialFileContent = harmonic_potential_1D())
        true
    end
end

"""
solver=lobpcg only ever supports non-periodic (real symmetric) problems, but
checkSolver runs before periodicity/k-points are known and so cannot reject
this combination itself - the check has to happen in setupSystem, once
system.reciprocal is set, and before main.jl's unconditional
`rm(files.eigenvalueFileName)`, so a user's existing results survive a
misconfigured re-run rather than being deleted ahead of the error.
"""
function test_lobpcg_periodic_rejected_before_side_effects()
    mktempdir() do tmp
        cd(tmp) do
            write("input.in", """
                               potential-file = potential.dat
                               solver         = lobpcg
                               periodic       = true
                               band-structure = on
                               k-points       = 10
                               mass-unit      = m_e
                               """)
            write("potential.dat", harmonic_potential_1D())
            write("eigenvalues.dat", "sentinel - must not be deleted by a rejected run")

            @test_throws ArgumentError @suppress Numerov.numerov("input.in")
            @test read("eigenvalues.dat", String) == "sentinel - must not be deleted by a rejected run"
        end
    end
end

function test_inputValidation()
    test_readInputFile_errors()
    test_readInputFile_valid()
    test_checkPotential_errors()
    test_checkPotential_accepted()
    test_checkSystem_errors()
    test_checkSystem_accepted()
    test_checkOutput()
    test_numerov_inputValidation()
    test_lobpcg_periodic_rejected_before_side_effects()
end
