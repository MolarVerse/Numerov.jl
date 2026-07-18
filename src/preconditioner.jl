"""
    KineticPreconditioner

Tensor-product approximation to the inverse of the shifted kinetic operator,
used to precondition the LOBPCG eigensolver.

The kinetic energy is approximated by the Kronecker sum of the per-dimension
1D operators `t_d = -Δ_d / (2 Δq_d²)`; its eigendecomposition factorizes into
the per-dimension eigenpairs, so `(T̃ + σI)⁻¹ x` is applied exactly with one
small dense eigenbasis transform per dimension - no factorization of the full
operator and therefore no fill-in.
"""
struct KineticPreconditioner
    Q     ::Vector{Matrix{Float64}}   # eigenbasis per dimension, d = 1..D
    denom ::Array{Float64}            # Σ_d λ_d + σ, shaped (n_D, ..., n_1)
end

"""
1D kinetic matrix `-Δ/(2 Δq²)` for one dimension, built with the same stencil
machinery as the full operator.
"""
function kinetic_1d(n::Int, spacing::Float64, periodic::Bool, stencil::Int)
    pot1              = Potential()
    pot1.dimension    = 1
    pot1.n_datapoints = [n]
    pot1.periodic     = [periodic]

    sys1              = System()
    sys1.n_datapoints = pot1.n_datapoints
    sys1.periodic     = pot1.periodic
    sys1.reciprocal   = false
    sys1.stencil      = stencil
    sys1.stencilΔ     = stencil
    sys1.stencil∇     = stencil

    buildΔ(sys1, pot1)

    return Symmetric(0.5 .* (-Matrix(sys1.Δ) ./ spacing^2))
end

"""
    KineticPreconditioner(potential, system; σ = 1.0)

Build the preconditioner for the given problem. `σ > 0` keeps the operator
safely positive definite; its exact value only affects convergence speed.
"""
function KineticPreconditioner(potential::Potential, system::System; σ::Float64 = 1.0)
    D    = potential.dimension
    dims = potential.n_datapoints

    Q       = Vector{Matrix{Float64}}(undef, D)
    lambdas = Vector{Vector{Float64}}(undef, D)
    for d in 1:D
        e = eigen(kinetic_1d(dims[d], potential.intervall[d], potential.periodic[d], system.stencilΔ))
        Q[d], lambdas[d] = e.vectors, e.values
    end

    # flattened grid ordering: the LAST dimension varies fastest, so array
    # axis k corresponds to dimension D + 1 - k
    denom = zeros(reverse(Tuple(dims)))
    for idx in CartesianIndices(denom)
        acc = σ
        for k in 1:D
            acc += lambdas[D + 1 - k][idx[k]]
        end
        denom[idx] = acc
    end

    return KineticPreconditioner(Q, denom)
end

"""
Apply `(T̃ + σI)⁻¹` to a flattened grid vector via per-dimension eigenbasis
transforms.
"""
function apply_preconditioner!(y::AbstractVector, P::KineticPreconditioner, x::AbstractVector)
    D    = length(P.Q)
    dims = size(P.denom)
    X    = reshape(copy(convert(Vector{Float64}, x)), dims)

    X = transform_modes(X, P, adjoint)
    X ./= P.denom
    X = transform_modes(X, P, identity)

    y .= vec(X)
    return y
end

"Multiply every mode of `X` by `op(Q_d)` for its dimension's eigenbasis."
function transform_modes(X::AbstractArray, P::KineticPreconditioner, op)
    D = length(P.Q)
    for k in 1:D
        d    = D + 1 - k                     # dimension of array axis k
        perm = (k, setdiff(1:D, k)...)
        Xp   = permutedims(X, perm)
        sz   = size(Xp)
        Xm   = op(P.Q[d]) * reshape(Xp, sz[1], :)
        X    = permutedims(reshape(Xm, sz), invperm(collect(perm)))
    end
    return X
end

LinearAlgebra.ldiv!(P::KineticPreconditioner, x::AbstractVector) =
    apply_preconditioner!(x, P, copy(x))

function LinearAlgebra.ldiv!(P::KineticPreconditioner, X::AbstractMatrix)
    for j in axes(X, 2)
        col = view(X, :, j)
        apply_preconditioner!(col, P, copy(col))
    end
    return X
end

LinearAlgebra.ldiv!(y::AbstractVector, P::KineticPreconditioner, x::AbstractVector) =
    apply_preconditioner!(y, P, x)

function LinearAlgebra.ldiv!(Y::AbstractMatrix, P::KineticPreconditioner, X::AbstractMatrix)
    for j in axes(X, 2)
        apply_preconditioner!(view(Y, :, j), P, view(X, :, j))
    end
    return Y
end
