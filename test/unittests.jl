using Numerov.MyUnits

include("unittests/test_1Dlaplace.jl")
include("unittests/test_1Dnabla.jl")
include("unittests/test_2Dlaplace.jl")
include("unittests/test_2Dnabla.jl")

include("unittests/test_k_paths.jl")

include("unittests/test_internalUnits.jl")

function unittests()
    @testset "test 1D laplace operator" test_1DΔ()
    @testset "test 1D nabla   operator" test_1D∇()
    @testset "test 2D laplace operator" test_2DΔ()
    @testset "test 2D nabla   operator" test_2D∇()

    k_intervalls = [1.0, 1.0, 1.0]
    n_kpoints    = 3

    @testset "test Γ_X" test_calc_Γ_X(k_intervalls, n_kpoints)
    @testset "test X_M" test_calc_X_M(k_intervalls, n_kpoints)
    @testset "test M_Γ" test_calc_M_Γ(k_intervalls, n_kpoints)
    @testset "test Γ_R" test_calc_Γ_R(k_intervalls, n_kpoints)
    @testset "test R_X" test_calc_R_X(k_intervalls, n_kpoints)
    @testset "test M_R" test_calc_M_R(k_intervalls, n_kpoints)

    @testset "test 3D k-path cube"            test_get_kpoints_3D_cube([1.0, 1.0, 1.0], n_kpoints)
    @testset "test 3D k-path noncube xysym"   test_get_kpoints_3D_noncube_xysym([1.0, 1.0, 2.0], n_kpoints)
    @testset "test 3D k-path noncube xzsym"   test_get_kpoints_3D_noncube_xzsym([1.0, 2.0, 1.0], n_kpoints)
    @testset "test 3D k-path noncube yzsym"   test_get_kpoints_3D_noncube_yzsym([2.0, 1.0, 1.0], n_kpoints)
    @testset "test 3D k-path noncube"         test_get_kpoints_3D_noncube([1.0, 2.0, 3.0], n_kpoints)
    @testset "test 3D k-path sym consistency" test_get_kpoints_3D_sym_consistency(n_kpoints)

    @testset "test internal units" test_internalUnits()
end

# function test_convert_to_internalUnits()
    
# end