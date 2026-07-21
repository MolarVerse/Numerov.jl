####################################################################
#                                                                  #
# tests for the programmatic API (src/api.jl)                      #
#                                                                  #
#   - equivalence with the file pipeline on the shipped test cases #
#   - unit handling (explicit units and Unitful quantities)        #
#   - single-k solves reproduce band-structure rows                #
#   - input validation error paths                                 #
#                                                                  #
####################################################################

api_testcase_path(parts...) = joinpath(@__DIR__, "..", "testsets", parts...)

api_maxabs(a, b) = maximum(abs.(a .- b))

"""
Copy all input files of `test/testsets/<case>/` into a fresh temporary
directory, run the file pipeline (`Numerov.numerov`) there and return the
parsed numeric output files (`eigenvalues.dat` and, if written,
`bandstructure.dat`).
"""
function api_run_file_pipeline(case::String)
    case_path = api_testcase_path(case)
    mktempdir() do tmp
        for file in readdir(case_path)
            cp(joinpath(case_path, file), joinpath(tmp, file))
        end
        cd(tmp) do
            @suppress Numerov.numerov("input.in")
            outputs = Dict{String, Matrix{Float64}}()
            for name in ("eigenvalues.dat", "bandstructure.dat")
                isfile(name) && (outputs[name] = Float64.(readdlm(name; comments = true)))
            end
            outputs
        end
    end
end

"""
Read a 1D `potential.dat` of a test case into `(coords, V)`.
"""
function api_read_potential_1D(case::String)
    data = readdlm(api_testcase_path(case, "potential.dat"); comments = true)
    return data[:, 1], data[:, 2]
end

function api_check_normalization(states, volume_element)
    for j in axes(states, 2)
        @test sum(abs2, view(states, :, j)) * volume_element ≈ 1.0 atol = 1.0e-10
    end
end

####################################################################
#                                                                  #
# 1. equivalence with the file pipeline                            #
#                                                                  #
# the file pipeline prints eigenvalues referenced to min(V) while  #
# the API returns absolute energies, so minimum(V) is subtracted   #
# from the API energies before comparing (== 0 for all three       #
# cases, kept for clarity)                                         #
#                                                                  #
####################################################################

function test_api_equivalence()

    ################################################################
    # 1D harmonic oscillator, non-periodic                         #
    # input.in: stencil 9, mass-unit me, coord-unit bohr,          #
    # n-eigenvalues 10, defaults: potential-unit hartree, mass 1.0 #
    ################################################################

    x, V = api_read_potential_1D("1DHarmonicOscillator")

    solution = solve_schrodinger(V, x;
        mass = 1.0, n_eigenvalues = 10, stencil = 9,
        potential_unit = u"hartree", coord_unit = u"bohr", mass_unit = u"m_e")

    reference = vec(api_run_file_pipeline("1DHarmonicOscillator")["eigenvalues.dat"])

    @test length(solution.energies) == 10
    @test issorted(solution.energies)
    @test api_maxabs(solution.energies .- minimum(V), reference) < 1.0e-7
    @test solution.kpoint === nothing
    @test eltype(solution.states) == Float64

    # states are column-normalized with the spacing in coordinate units (bohr)
    api_check_normalization(solution.states, x[2] - x[1])

    ################################################################
    # 2D water potential energy surface, non-periodic              #
    # input.in: stencil 9, mass-unit unit (-> u), coord-unit       #
    # angstrom, n-eigenvalues 5, datapoints 66 53, potential-unit  #
    # kcal/mol, default mass 1.0                                   #
    ################################################################

    data = readdlm(api_testcase_path("2DWater", "potential.dat"); comments = true)

    x2 = unique(data[:, 1])
    y2 = unique(data[:, 2])

    # the file rows run x-outer / y-inner ...
    @test data[:, 1] == repeat(x2, inner = length(y2))
    @test data[:, 2] == repeat(y2, outer = length(x2))

    # ... so a column-major reshape gives (iy, ix) and a permutedims V[ix, iy]
    V2 = permutedims(reshape(data[:, 3], length(y2), length(x2)))
    @test size(V2) == (66, 53) # "datapoints = 66 53" in input.in

    solution2 = solve_schrodinger(V2, (x2, y2);
        mass = 1.0, n_eigenvalues = 5, stencil = 9,
        potential_unit = u"kcalpermol", coord_unit = u"angstrom", mass_unit = u"u")

    reference2 = vec(api_run_file_pipeline("2DWater")["eigenvalues.dat"])

    @test length(solution2.energies) == 5
    @test api_maxabs(solution2.energies .- minimum(V2), reference2) < 1.0e-7

    api_check_normalization(solution2.states, (x2[2] - x2[1]) * (y2[2] - y2[1]))

    ################################################################
    # 1D Kronig-Penney band structure, periodic                    #
    # input.in: stencil 9, reduced-mass 1.0, mass-unit me,         #
    # coord-unit angstrom, n-eigenvalues 5, potential-unit ev,     #
    # periodic true, band-structure on, k-points 10                #
    ################################################################

    xkp, Vkp = api_read_potential_1D("1DKronigPenney")

    bands = band_structure(Vkp, xkp;
        n_kpoints = 10, mass = 1.0, periodic = true, n_eigenvalues = 5,
        stencil = 9, potential_unit = u"eV", coord_unit = u"angstrom", mass_unit = u"m_e")

    reference_bs = api_run_file_pipeline("1DKronigPenney")["bandstructure.dat"]

    @test size(reference_bs) == (10, 6)
    @test size(bands.energies) == (10, 5)
    @test length(bands.kpoints) == 10
    @test bands.kpoints[1] == [0.0] # the path starts at Γ
    @test bands.kpath[1] == 0.0
    @test issorted(bands.kpath)

    # the distance column of bandstructure.dat is in internal mass-weighted
    # units (bohr⁻¹ for mass = 1 mₑ) while the API kpath is in coordinate
    # units (Å⁻¹), so only the length conversion is needed here
    kpath_internal = ustrip.(u"bohr^-1", bands.kpath .* u"angstrom^-1")
    @test api_maxabs(kpath_internal, reference_bs[:, 1]) < 1.0e-6

    @test api_maxabs(bands.energies .- minimum(Vkp), reference_bs[:, 2:end]) < 1.0e-6
end

####################################################################
#                                                                  #
# 2. unit handling                                                 #
#                                                                  #
####################################################################

function test_api_units()
    x_bohr    = collect(range(-10.0, 10.0; length = 201))
    V_hartree = 0.5 .* x_bohr .^ 2

    reference = solve_schrodinger(V_hartree, x_bohr; n_eigenvalues = 6)

    # analytic sanity check: E_n = n + 1/2 for ω = m = ħ = 1
    @test api_maxabs(reference.energies, collect(0:5) .+ 0.5) < 1.0e-5

    # the same problem from pre-converted plain arrays with explicit units
    V_ev  = ustrip.(u"eV", V_hartree .* u"hartree")
    x_ang = ustrip.(u"angstrom", x_bohr .* u"bohr")

    solution_ev = solve_schrodinger(V_ev, x_ang; n_eigenvalues = 6,
        potential_unit = u"eV", coord_unit = u"angstrom")

    @test api_maxabs(ustrip.(u"hartree", solution_ev.energies .* u"eV"),
                     reference.energies) < 1.0e-8

    # states are normalized in the coordinate unit actually used (angstrom)
    api_check_normalization(solution_ev.states, x_ang[2] - x_ang[1])

    # Unitful convenience method: quantities in, hartree/bohr out
    solution_quantity = solve_schrodinger(V_ev .* u"eV", x_ang .* u"angstrom";
        n_eigenvalues = 6)

    @test api_maxabs(solution_quantity.energies, reference.energies) < 1.0e-8

    api_check_normalization(solution_quantity.states, x_bohr[2] - x_bohr[1])
end

####################################################################
#                                                                  #
# 3. a single-k solve reproduces the band-structure row            #
#                                                                  #
####################################################################

function test_api_single_k()
    n_points = 120
    dx       = 0.05
    L        = n_points * dx

    x = collect(range(0.0; step = dx, length = n_points))
    V = 0.3 .* (1 .- cos.(2π .* x ./ L))

    bands = band_structure(V, x; n_kpoints = 6, n_eigenvalues = 4, periodic = true)

    @test length(bands.kpoints) == 6
    @test issorted(bands.kpath)

    index    = 4
    solution = solve_schrodinger(V, x;
        periodic = true, k = bands.kpoints[index], n_eigenvalues = 4)

    @test solution.kpoint == bands.kpoints[index]
    @test eltype(solution.states) <: Complex
    @test api_maxabs(solution.energies, bands.energies[index, :]) < 1.0e-9
end

####################################################################
#                                                                  #
# 4. error paths of the input validation                           #
#                                                                  #
####################################################################

function test_api_errors()
    x  = collect(range(0.0, 1.0; length = 32))
    V1 = x .^ 2
    V2 = [xi^2 + yi^2 for xi in x, yi in x]

    z  = x[1:12]
    V3 = [xi + yi + zi for xi in z, yi in z, zi in z]

    # wrong number of coordinate axes
    @test_throws ArgumentError solve_schrodinger(V2, (x,))
    @test_throws ArgumentError solve_schrodinger(V1, (x, x))

    # more than 3 dimensions
    @test_throws ArgumentError solve_schrodinger(ones(2, 2, 2, 2), (1:2, 1:2, 1:2, 1:2))

    # axis length does not match the potential
    @test_throws ArgumentError solve_schrodinger(V1, x[1:31])

    # descending axis
    @test_throws ArgumentError solve_schrodinger(V1, reverse(x))

    # unequally spaced axis (still ascending)
    x_uneven     = copy(x)
    x_uneven[5] += 0.01
    @test_throws ArgumentError solve_schrodinger(V1, x_uneven)

    # equal plain spacings but different masses -> the mass-weighted
    # spacings dx * sqrt(m) differ between the dimensions
    @test_throws ArgumentError solve_schrodinger(V2, (x, x); mass = [1.0, 4.0])

    # a k-point requires at least one periodic dimension
    @test_throws ArgumentError solve_schrodinger(V1, x; k = [0.1])

    # unknown solver name
    @test_throws ArgumentError solve_schrodinger(V1, x; solver = :cuda)

    # n_eigenvalues has to be positive
    @test_throws ArgumentError solve_schrodinger(V1, x; n_eigenvalues = 0)

    # even stencil sizes do not exist
    @test_throws ArgumentError solve_schrodinger(V1, x; stencil = 4)

    # the nabla stencil has no 13-point variant (checked also on 3D input)
    @test_throws ArgumentError solve_schrodinger(V3, (z, z, z); stencil_nabla = 13)

    # n_eigenvalues too close to the grid size is rejected up front for every
    # solver that shares arpack's nev < N requirement - :arpack directly,
    # :lobpcg via its arpack fallback, and :krylov - not just :arpack; before
    # this guard, :lobpcg raised an opaque BoundsError deep in solve() instead
    for solver in (:arpack, :lobpcg, :krylov)
        @test_throws ArgumentError solve_schrodinger(V1, x; n_eigenvalues = length(x), solver = solver)
    end
end
