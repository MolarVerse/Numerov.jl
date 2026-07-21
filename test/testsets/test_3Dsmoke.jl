"""
    generate_3DHO_files(dir, n_points, xmin, xmax, stencil)

Write the potential and input file of a 3D isotropic harmonic oscillator
V(x,y,z) = 0.5⋅(x² + y² + z²) (atomic units, ω = 1) on a cubic grid with
`n_points` points per dimension spanning [`xmin`, `xmax`] into `dir`.

The grid-point ordering matches what `readPotential` expects for 3D
potentials (x outermost loop, z innermost - "x y z V" per line, compare
`examples/3DKronigPenney/potential.dat`).
"""
function generate_3DHO_files(dir::String, n_points::Int, xmin::Float64, xmax::Float64, stencil::Int)

    coords = range(xmin, xmax; length=n_points)

    open(joinpath(dir, "potential.dat"), "w") do file
        for x in coords, y in coords, z in coords
            println(file, x, "   ", y, "   ", z, "   ", 0.5 * (x^2 + y^2 + z^2))
        end
    end

    open(joinpath(dir, "input.in"), "w") do file
        println(file, "potential-file = potential.dat")
        println(file, "potential-unit = hartree")
        println(file, "coord-unit = bohr")
        println(file, "mass-unit = me")
        println(file, "reduced-mass = 1.0")
        println(file, "datapoints = $n_points, $n_points, $n_points")
        println(file, "n-eigenvalues = 4")
        println(file, "stencil = $stencil")
    end

end

"""
    test_3Dsmoke()

Analytic smoke test for the 3D Laplacian stencils: solves a small 3D
isotropic harmonic oscillator for every supported stencil size and checks
the eigenvalues against the analytic spectrum Eₙ = (n_x + n_y + n_z + 3/2) ħω,
i.e. E₀ = 1.5 Eₕ and the threefold degenerate E₁ = 2.5 Eₕ. Also checks that
the unsupported 3- and 13-point stencils throw an `ArgumentError` in 3D.
"""
function test_3Dsmoke()

    n_points = 15
    xmin     = -5.5
    xmax     =  5.5

    ##################################################################################
    #                                                                                #
    # tolerances calibrated empirically on the 15³ grid with roughly a 2x safety     #
    # margin over the observed discretization error per stencil size:               #
    #                                                                                #
    #   stencil => (atol E₀, atol E₁)   observed errors: (E₀, E₁)                   #
    #        5  => 0.0055, 0.017                                                     #
    #        7  => 0.0024, 0.0083                                                    #
    #        9  => 0.0015, 0.0058                                                    #
    #       11  => 0.0011, 0.0041                                                    #
    #                                                                                #
    ##################################################################################

    tolerances = Dict(
         5 => (0.011 , 0.034 ),
         7 => (0.005 , 0.017 ),
         9 => (0.0032, 0.012 ),
        11 => (0.0022, 0.0085),
    )

    for stencil in (5, 7, 9, 11)
        mktempdir() do tmp

            generate_3DHO_files(tmp, n_points, xmin, xmax, stencil)

            cd(tmp) do
                @suppress Numerov.numerov("input.in")

                eigenvalues = vec(readdlm("eigenvalues.dat"; comments=true))

                atol_E0, atol_E1 = tolerances[stencil]

                @test length(eigenvalues) == 4
                @test issorted(eigenvalues)

                @test eigenvalues[1] ≈ 1.5 atol = atol_E0

                for i in 2:4
                    @test eigenvalues[i] ≈ 2.5 atol = atol_E1
                end

                # first excited level has to be threefold degenerate on the cubic grid
                @test eigenvalues[4] - eigenvalues[2] < 1.0e-6
            end

        end
    end

    #####################################################################
    #                                                                   #
    # 3- and 13-point stencils are not implemented for 3d calculations #
    #                                                                   #
    #####################################################################

    # the lobpcg solver must reproduce the arpack analytic result - and must
    # actually do so via lobpcg itself, not via solveWrapper's automatic
    # retry/escalation-to-arpack fallback, which would make `r.energies` look
    # correct even if lobpcg silently failed on this genuinely (3-fold)
    # degenerate cluster and Arpack quietly rescued it. Seed the RNG for
    # reproducibility and assert no fallback/escalation warning fired, so
    # this test actually proves lobpcg converged rather than merely that the
    # pipeline's overall answer is correct.
    let x = range(xmin, xmax; length = n_points)
        V = [0.5 * (a^2 + b^2 + c^2) for a in x, b in x, c in x]
        (a0, a1) = (0.0032, 0.012)

        Random.seed!(1)
        local r
        logs, _ = Test.collect_test_logs() do
            r = solve_schrodinger(V, (x, x, x); n_eigenvalues = 4, solver = :lobpcg)
        end
        rescued = any(
            l -> occursin("falling back to arpack", l.message) ||
                 occursin("exceed the residual tolerance", l.message),
            logs,
        )
        @test !rescued
        @test isapprox(r.energies[1], 1.5; atol = a0)
        @test all(isapprox.(r.energies[2:4], 2.5; atol = a1))
    end

    for stencil in (3, 13)
        mktempdir() do tmp

            generate_3DHO_files(tmp, n_points, xmin, xmax, stencil)

            cd(tmp) do
                @test_throws ArgumentError @suppress Numerov.numerov("input.in")
            end

        end
    end

end
