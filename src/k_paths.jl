function get_kpoints_1D(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    return map(x -> Tuple(x[1]), calc_Γ_X(k_intervalls, n_kpoints-1))
end

function get_kpoints_2D(k_intervalls::Vector{Float64}, n_kpoints::Int64)
    if isapprox(k_intervalls[1], k_intervalls[2], atol=1e-9)
        return get_kpoints_2D_square(k_intervalls, n_kpoints)
    else
        return get_kpoints_2D_rectangle(k_intervalls, n_kpoints)
    end

end

function get_kpoints_2D_square(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    X_S = calc_X_S(k_intervalls, n_kpoints-1)
    
    S_Γ = calc_S_Γ(k_intervalls, n_kpoints-1)

    minimal_kpath = rle(vcat(Γ_X, X_S, S_Γ))[1]

    return map(x -> Tuple(x[1:2]), minimal_kpath)
end

function get_kpoints_2D_rectangle(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    X_S = calc_X_S(k_intervalls, n_kpoints-1)
    
    S_Y = calc_S_Y(k_intervalls, n_kpoints-1)

    Y_Γ = calc_Y_Γ(k_intervalls, n_kpoints-1)

    minimal_kpath = rle(vcat(Γ_X, X_S, S_Y, Y_Γ))[1]

    return map(x -> Tuple(x[1:2]), minimal_kpath)
end

function get_kpoints_3D(k_intervalls, n_kpoints)
    if isapprox(k_intervalls[1], k_intervalls[2], atol=1e-9)

        if isapprox(k_intervalls[1], k_intervalls[3], atol=1e-9)

            get_kpoints_3D_cube(k_intervalls, n_kpoints)

        else

            get_kpoints_3D_noncube_xysym(k_intervalls, n_kpoints)

        end

    elseif isapprox(k_intervalls[1], k_intervalls[3], atol=1e-9)

        get_kpoints_3D_noncube_xzsym(k_intervalls, n_kpoints)

    elseif isapprox(k_intervalls[2], k_intervalls[3], atol=1e-9)

        get_kpoints_3D_noncube_yzsym(k_intervalls, n_kpoints)

    else

        get_kpoints_3D_noncube(k_intervalls, n_kpoints)

    end
end

function get_kpoints_3D_cube(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    X_S = calc_X_S(k_intervalls, n_kpoints-1)
    
    S_Γ = calc_S_Γ(k_intervalls, n_kpoints-1)

    Γ_R = calc_Γ_R(k_intervalls, n_kpoints-1)

    R_X = calc_R_X(k_intervalls, n_kpoints-1)
    
    S_R = calc_S_R(k_intervalls, n_kpoints-1)
    
    return rle(vcat(Γ_X, X_S, S_Γ, Γ_R, R_X, S_R))[1]
end

function get_kpoints_3D_noncube_xysym(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    X_S = calc_X_S(k_intervalls, n_kpoints-1)
    
    S_Γ = calc_S_Γ(k_intervalls, n_kpoints-1)

    Γ_Z = calc_Γ_Z(k_intervalls, n_kpoints-1)

    Z_T = calc_Z_T(k_intervalls, n_kpoints-1)

    T_R = calc_T_R(k_intervalls, n_kpoints-1)

    Z_R = calc_Z_R(k_intervalls, n_kpoints-1)

    X_T = calc_X_T(k_intervalls, n_kpoints-1)

    S_R = calc_S_R(k_intervalls, n_kpoints-1)
    
    return rle(vcat(Γ_X, X_S, S_Γ, Γ_Z, Z_T, T_R, Z_R, T_R, Z_R, X_T, S_R))[1]
end

function get_kpoints_3D_noncube_xzsym(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    X_T = calc_X_T(k_intervalls, n_kpoints-1)
    
    T_Γ = calc_T_Γ(k_intervalls, n_kpoints-1)

    Γ_Y = calc_Γ_Y(k_intervalls, n_kpoints-1)

    Y_S = calc_Y_S(k_intervalls, n_kpoints-1)

    S_R = calc_S_R(k_intervalls, n_kpoints-1)

    Y_R = calc_Y_R(k_intervalls, n_kpoints-1)

    X_S = calc_X_S(k_intervalls, n_kpoints-1)

    T_R = calc_T_R(k_intervalls, n_kpoints-1)
    
    return rle(vcat(Γ_X, X_T, T_Γ, Γ_Y, Y_S, S_R, Y_R, T_R, Z_R, X_S, S_R, T_R))[1]
end

function get_kpoints_3D_noncube_yzsym(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_Y = calc_Γ_Y(k_intervalls, n_kpoints-1)

    Y_U = calc_Y_U(k_intervalls, n_kpoints-1)
    
    U_Γ = calc_U_Γ(k_intervalls, n_kpoints-1)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    Z_S = calc_X_S(k_intervalls, n_kpoints-1)

    S_R = calc_S_R(k_intervalls, n_kpoints-1)

    X_R = calc_X_R(k_intervalls, n_kpoints-1)

    Y_S = calc_Y_S(k_intervalls, n_kpoints-1)

    U_R = calc_U_R(k_intervalls, n_kpoints-1)
    
    return rle(vcat(Γ_Y, Y_U, U_Γ, Γ_X, X_S, S_R, X_R, Y_S, U_R))[1]
end

function get_kpoints_3D_noncube(k_intervalls, n_kpoints)
    
    Γ_X = calc_Γ_X(k_intervalls, n_kpoints-1)

    X_S = calc_X_S(k_intervalls, n_kpoints-1)

    S_Y = calc_S_Y(k_intervalls, n_kpoints-1)

    Y_Γ = calc_Y_Γ(k_intervalls, n_kpoints-1)

    Γ_Z = calc_Γ_Z(k_intervalls, n_kpoints-1)

    Z_U = calc_Z_U(k_intervalls, n_kpoints-1)

    U_R = calc_U_R(k_intervalls, n_kpoints-1)

    R_T = calc_R_T(k_intervalls, n_kpoints-1)

    T_Z = calc_T_Z(k_intervalls, n_kpoints-1)

    Y_T = calc_Y_T(k_intervalls, n_kpoints-1)

    U_X = calc_U_X(k_intervalls, n_kpoints-1)

    S_R = calc_S_R(k_intervalls, n_kpoints-1)

    return rle(vcat(Γ_X, X_S, S_Y, Y_Γ, Γ_Z, Z_U, U_R, R_T, T_Z, Y_T, U_X, S_R))[1]
end

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

calc_Γ_Z(k, n) = [(   0.0,    0.0, i*k[3]) for i in 0:n]

calc_Γ_Y(k, n) = [(   0.0, i*k[3],    0.0) for i in 0:n]

calc_Z_U(k, n) = [(   0.0, i*k[2], n*k[3]) for i in 0:n]

calc_Y_U(k, n) = [(   0.0, n*k[2], i*k[3]) for i in 0:n]

calc_Γ_X(k, n) = [(i*k[1],    0.0,    0.0) for i in 0:n]

calc_Z_T(k, n) = [(i*k[1],    0.0, n*k[3]) for i in 0:n]

calc_Γ_R(k, n) = [(i*k[1], i*k[2], i*k[3]) for i in 0:n]

calc_Y_S(k, n) = [(i*k[1], n*k[2],    0.0) for i in 0:n]

calc_Y_R(k, n) = [(i*k[1], n*k[2], i*k[2]) for i in 0:n]

calc_U_R(k, n) = [(i*k[1], n*k[2], n*k[2]) for i in 0:n]

calc_X_T(k, n) = [(n*k[1],    0.0, i*k[3]) for i in 0:n]

calc_X_S(k, n) = [(n*k[1], i*k[2],    0.0) for i in 0:n]

calc_X_R(k, n) = [(n*k[1], i*k[2], i*k[3]) for i in 0:n]

calc_T_R(k, n) = [(n*k[1], i*k[2], n*k[3]) for i in 0:n]

calc_S_R(k, n) = [(n*k[1], n*k[2], i*k[3]) for i in 0:n]

###########################################################

calc_Y_Γ(k, n) = [(   0.0, i*k[2],    0.0) for i in n:-1:0]

calc_U_Γ(k, n) = [(   0.0, i*k[2], i*k[3]) for i in n:-1:0]

calc_T_Γ(k, n) = [(i*k[1],    0.0, i*k[3]) for i in n:-1:0]

calc_T_Z(k, n) = [(i*k[1],    0.0, n*k[3]) for i in n:-1:0]

calc_S_Γ(k, n) = [(i*k[1], i*k[2],    0.0) for i in n:-1:0]

calc_Z_R(k, n) = [(i*k[1], i*k[2], n*k[3]) for i in n:-1:0]

calc_S_Y(k, n) = [(i*k[1], n*k[2],    0.0) for i in n:-1:0]

calc_R_X(k, n) = [(n*k[1], i*k[2], i*k[3]) for i in n:-1:0]

calc_R_T(k, n) = [(n*k[1], i*k[2], n*k[3]) for i in n:-1:0]

###########################################################

calc_Y_T(k, n) = [(i*k[1], (n-i)*k[2],     i*k[3]) for i in 0:n]

calc_U_X(k, n) = [(i*k[1], (n-i)*k[2], (n-i)*k[3]) for i in 0:n]