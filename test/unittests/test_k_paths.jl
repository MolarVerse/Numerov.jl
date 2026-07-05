        ######################################
       #.                                   ##
      # .                                  # #
     #  .                                 #  #
    #   .              Z ----------------T   #
   #    .             --                #-   # 
  #     .            - -               # -   # 
 #      .           -  -              #  -   # 
################## U ############### R   -   # 
#       .          -   -             #   -   # 
#       .          -   -             #   -   # 
#       .          -   Γ-------------#---X   # 
#       .          -  -              #  -    # 
#       .          - -               # -     # 
#       ...........--................#-......# 
#      .           Y-----------------S      #  
#     .                              #     #   
#    .                               #    #    
#   .                                #   #     
#  .                                 #  #  
# .                                  # #   
#.                                   ##                     
###################################### 

function test_calc_Γ_X(k_intervalls, n_kpoints)
    @test Numerov.calc_Γ_X(k_intervalls, n_kpoints-1) == [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
end

function test_calc_X_M(k_intervalls, n_kpoints)
    @test Numerov.calc_X_S(k_intervalls, n_kpoints-1) == [(2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (2.0, 2.0, 0.0)]
end

function test_calc_M_Γ(k_intervalls, n_kpoints)
    @test Numerov.calc_S_Γ(k_intervalls, n_kpoints-1) == [(2.0, 2.0, 0.0), (1.0, 1.0, 0.0), (0.0, 0.0, 0.0)]
end

function test_calc_Γ_R(k_intervalls, n_kpoints)
    @test Numerov.calc_Γ_R(k_intervalls, n_kpoints-1) == [(0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2.0, 2.0, 2.0)]
end

function test_calc_R_X(k_intervalls, n_kpoints)
    @test Numerov.calc_R_X(k_intervalls, n_kpoints-1) == [(2.0, 2.0, 2.0), (2.0, 1.0, 1.0), (2.0, 0.0, 0.0)]
end

function test_calc_M_R(k_intervalls, n_kpoints)
    @test Numerov.calc_S_R(k_intervalls, n_kpoints-1) == [(2.0, 2.0, 0.0), (2.0, 2.0, 1.0), (2.0, 2.0, 2.0)]
end

function check_no_consecutive_duplicates(kpath)
    @test all(kpath[i] != kpath[i+1] for i in 1:length(kpath)-1)
end

function test_get_kpoints_3D_cube(k_intervalls, n_kpoints)
    kpath = Numerov.get_kpoints_3D(k_intervalls, n_kpoints)

    # Γ-X-S-Γ-R-X | S-R with Γ = (0,0,0), X = (2,0,0), S = (2,2,0), R = (2,2,2)
    @test kpath == [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                    (2.0, 1.0, 0.0), (2.0, 2.0, 0.0), (1.0, 1.0, 0.0),
                    (0.0, 0.0, 0.0), (1.0, 1.0, 1.0), (2.0, 2.0, 2.0),
                    (2.0, 1.0, 1.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0),
                    (2.0, 2.0, 1.0), (2.0, 2.0, 2.0)]

    check_no_consecutive_duplicates(kpath)
end

function test_get_kpoints_3D_noncube_xysym(k_intervalls, n_kpoints)
    kpath = Numerov.get_kpoints_3D(k_intervalls, n_kpoints)

    # Γ-X-S-Γ-Z-T-R-Z | X-T | S-R with X = (2,0,0), S = (2,2,0), Z = (0,0,4), T = (2,0,4), R = (2,2,4)
    @test kpath == [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                    (2.0, 1.0, 0.0), (2.0, 2.0, 0.0), (1.0, 1.0, 0.0),
                    (0.0, 0.0, 0.0), (0.0, 0.0, 2.0), (0.0, 0.0, 4.0),
                    (1.0, 0.0, 4.0), (2.0, 0.0, 4.0), (2.0, 1.0, 4.0),
                    (2.0, 2.0, 4.0), (1.0, 1.0, 4.0), (0.0, 0.0, 4.0),
                    (2.0, 0.0, 0.0), (2.0, 0.0, 2.0), (2.0, 0.0, 4.0),
                    (2.0, 2.0, 0.0), (2.0, 2.0, 2.0), (2.0, 2.0, 4.0)]

    check_no_consecutive_duplicates(kpath)
end

function test_get_kpoints_3D_noncube_xzsym(k_intervalls, n_kpoints)
    kpath = Numerov.get_kpoints_3D(k_intervalls, n_kpoints)

    # Γ-X-T-Γ-Y-S-R | Y-R | X-S | T-R with X = (2,0,0), T = (2,0,2), Y = (0,4,0), S = (2,4,0), R = (2,4,2)
    @test kpath == [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                    (2.0, 0.0, 1.0), (2.0, 0.0, 2.0), (1.0, 0.0, 1.0),
                    (0.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 4.0, 0.0),
                    (1.0, 4.0, 0.0), (2.0, 4.0, 0.0), (2.0, 4.0, 1.0),
                    (2.0, 4.0, 2.0), (0.0, 4.0, 0.0), (1.0, 4.0, 1.0),
                    (2.0, 4.0, 2.0), (2.0, 0.0, 0.0), (2.0, 2.0, 0.0),
                    (2.0, 4.0, 0.0), (2.0, 0.0, 2.0), (2.0, 2.0, 2.0),
                    (2.0, 4.0, 2.0)]

    check_no_consecutive_duplicates(kpath)
end

function test_get_kpoints_3D_noncube_yzsym(k_intervalls, n_kpoints)
    kpath = Numerov.get_kpoints_3D(k_intervalls, n_kpoints)

    # Γ-Y-U-Γ-X-S-R | X-R | Y-S | U-R with Y = (0,2,0), U = (0,2,2), X = (4,0,0), S = (4,2,0), R = (4,2,2)
    @test kpath == [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 2.0, 0.0),
                    (0.0, 2.0, 1.0), (0.0, 2.0, 2.0), (0.0, 1.0, 1.0),
                    (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (4.0, 0.0, 0.0),
                    (4.0, 1.0, 0.0), (4.0, 2.0, 0.0), (4.0, 2.0, 1.0),
                    (4.0, 2.0, 2.0), (4.0, 0.0, 0.0), (4.0, 1.0, 1.0),
                    (4.0, 2.0, 2.0), (0.0, 2.0, 0.0), (2.0, 2.0, 0.0),
                    (4.0, 2.0, 0.0), (0.0, 2.0, 2.0), (2.0, 2.0, 2.0),
                    (4.0, 2.0, 2.0)]

    check_no_consecutive_duplicates(kpath)
end

function test_get_kpoints_3D_noncube(k_intervalls, n_kpoints)
    kpath = Numerov.get_kpoints_3D(k_intervalls, n_kpoints)

    # Γ-X-S-Y-Γ-Z-U-R-T-Z | Y-T | U-X | S-R with X = (2,0,0), Y = (0,4,0), Z = (0,0,6),
    # S = (2,4,0), U = (0,4,6), T = (2,0,6), R = (2,4,6)
    @test kpath == [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0),
                    (2.0, 2.0, 0.0), (2.0, 4.0, 0.0), (1.0, 4.0, 0.0),
                    (0.0, 4.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 0.0),
                    (0.0, 0.0, 3.0), (0.0, 0.0, 6.0), (0.0, 2.0, 6.0),
                    (0.0, 4.0, 6.0), (1.0, 4.0, 6.0), (2.0, 4.0, 6.0),
                    (2.0, 2.0, 6.0), (2.0, 0.0, 6.0), (1.0, 0.0, 6.0),
                    (0.0, 0.0, 6.0), (0.0, 4.0, 0.0), (1.0, 2.0, 3.0),
                    (2.0, 0.0, 6.0), (0.0, 4.0, 6.0), (1.0, 2.0, 3.0),
                    (2.0, 0.0, 0.0), (2.0, 4.0, 0.0), (2.0, 4.0, 3.0),
                    (2.0, 4.0, 6.0)]

    check_no_consecutive_duplicates(kpath)
end

function test_get_kpoints_3D_sym_consistency(n_kpoints)
    yzsym = Numerov.get_kpoints_3D([2.0, 1.0, 1.0], n_kpoints)
    xzsym = Numerov.get_kpoints_3D([1.0, 2.0, 1.0], n_kpoints)
    xysym = Numerov.get_kpoints_3D([1.0, 1.0, 2.0], n_kpoints)

    # swapping the x and y axes maps the yzsym path onto the xzsym path point by point
    @test map(p -> (p[2], p[1], p[3]), yzsym) == xzsym

    # the cyclic axis permutation (x,y,z) -> (y,z,x) maps the k-points visited
    # by the yzsym path onto those visited by the xysym path
    @test Set(map(p -> (p[2], p[3], p[1]), yzsym)) == Set(xysym)
end