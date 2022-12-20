function get_kpoints_1D(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    return map(x -> Tuple(x[1]), calc_Γ_X(k_intervalls, n_kpoints))
end

function get_kpoints_2D(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints)

    X_M = calc_X_M(k_intervalls, n_kpoints)
    
    M_Γ = calc_M_Γ(k_intervalls, n_kpoints)

    minimal_kpath = rle(vcat(Γ_X, X_M, M_Γ))[1]

    return map(x -> Tuple(x[1:2]), minimal_kpath)
end

function get_kpoints_3D(k_intervalls::Vector{Float64}, n_kpoints::Int64)

    Γ_X = calc_Γ_X(k_intervalls, n_kpoints)

    X_M = calc_X_M(k_intervalls, n_kpoints)
    
    M_Γ = calc_M_Γ(k_intervalls, n_kpoints)

    Γ_R = calc_Γ_R(k_intervalls, n_kpoints)

    R_X = calc_R_X(k_intervalls, n_kpoints)
    
    M_R = calc_M_R(k_intervalls, n_kpoints)
    
    return rle(vcat(Γ_X, X_M, M_Γ, Γ_R, R_X, M_R))[1]
end

calc_Γ_X(k_intervalls, n_kpoints) = [(k_intervalls[1]*i, 0.0, 0.0) for i in 0:n_kpoints-1]

calc_X_M(k_intervalls, n_kpoints) = [((n_kpoints-1)*k_intervalls[1], i*k_intervalls[2], 0.0) for i in 0:n_kpoints-1]
    
calc_M_Γ(k_intervalls, n_kpoints) = [(i*k_intervalls[1], i*k_intervalls[2], 0.0) for i in n_kpoints-1:-1:0]

calc_Γ_R(k_intervalls, n_kpoints) = [(i*k_intervalls[1], i*k_intervalls[2], i*k_intervalls[3]) for i in 0:n_kpoints-1]

calc_R_X(k_intervalls, n_kpoints) = [((n_kpoints-1)*k_intervalls[1], i*k_intervalls[2], i*k_intervalls[3]) for i in n_kpoints-1:-1:0]
    
calc_M_R(k_intervalls, n_kpoints) = [((n_kpoints-1)*k_intervalls[1], (n_kpoints-1)*k_intervalls[2], k_intervalls[3]*i) for i in 0:n_kpoints-1]