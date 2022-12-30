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