####################################################################
#                                                                  #
# programmatic API: solve directly from arrays, no input files, no #
# output files, no global state                                    #
#                                                                  #
####################################################################

"""
    SchrodingerSolution

Result of [`solve_schrodinger`](@ref).

# Fields

- `energies::Vector{Float64}`: the lowest eigenvalues in the unit the
  potential was given in (absolute, i.e. referenced to the same zero as the
  input potential).
- `states::Matrix`: the corresponding eigenvectors as columns, in the same
  grid ordering as the flattened input potential array. Normalized so that
  `sum(abs2, state) * prod(spacings) == 1` with the spacings expressed in the
  unit the coordinates were given in. Real for non-periodic calculations,
  complex for periodic ones.
- `kpoint::Union{Nothing, Vector{Float64}}`: the k-point the equation was
  solved at (in inverse coordinate units), or `nothing` for non-periodic
  calculations.
"""
struct SchrodingerSolution{T <: Number}
    energies  ::Vector{Float64}
    states    ::Matrix{T}
    kpoint    ::Union{Nothing, Vector{Float64}}
end

"""
    BandStructure

Result of [`band_structure`](@ref).

# Fields

- `kpoints::Vector{Vector{Float64}}`: the k-points along the path through the
  high-symmetry points of the Brillouin zone, in inverse coordinate units.
- `kpath::Vector{Float64}`: cumulative distance along the k-path (inverse
  coordinate units), useful as the x axis of a band-structure plot.
- `energies::Matrix{Float64}`: band energies, one row per k-point and one
  column per band, in the unit the potential was given in (absolute).
"""
struct BandStructure
    kpoints ::Vector{Vector{Float64}}
    kpath   ::Vector{Float64}
    energies::Matrix{Float64}
end

# accepted solver names for the keyword interface
const SOLVER_NAMES = Dict(
    :arpack => ARPACK,
    :krylov => KRYLOV,
    :lu     => LU,
)

"""
    solve_schrodinger(V, coords; kwargs...) -> SchrodingerSolution

Solve the time-independent Schrödinger equation for the potential `V` given
on a regular grid and return the lowest eigenpairs, without reading or
writing any files.

# Arguments

- `V::AbstractArray{<:Real}`: the potential on the grid, as a 1-, 2- or
  3-dimensional array. `V[i, j, k]` is the value at
  `(coords[1][i], coords[2][j], coords[3][k])`.
- `coords`: the grid axes — one `AbstractVector`/`AbstractRange` per
  dimension of `V` (a single vector may be passed directly for 1D). Each axis
  must be ascending and equally spaced.

# Keyword arguments

- `mass = 1.0`: reduced mass, a scalar or one value per dimension.
- `periodic = false`: periodic boundary conditions, a scalar or one flag per
  dimension.
- `n_eigenvalues = 5`: number of eigenpairs to compute.
- `k = nothing`: solve at this k-point (a vector in inverse coordinate
  units); requires at least one periodic dimension. `nothing` solves the
  real, non-periodic problem.
- `stencil = 9`: finite-difference stencil size (`3`, `5`, `7`, `9`, `11` or
  `13`); `stencil_laplace` and `stencil_nabla` override it individually.
- `solver = :arpack`: eigensolver backend (`:arpack`, `:krylov` or `:lu`).
- `potential_unit = u"hartree"`, `coord_unit = u"bohr"`, `mass_unit = u"m_e"`:
  units of the inputs; energies are returned in `potential_unit`.

The grid spacings must be equal across dimensions after mass weighting, i.e.
`step(coords[i]) * sqrt(mass[i])` must be the same for every dimension — this
is a restriction of the underlying discretization.

# Example

```julia
x = range(-5.0, 5.0; length = 201)
result = solve_schrodinger(0.5 .* x .^ 2, x; n_eigenvalues = 3)
result.energies  # ≈ [0.5, 1.5, 2.5]
```
"""
function solve_schrodinger(V::AbstractArray{<:Real}, coords;
                           mass = 1.0,
                           periodic = false,
                           n_eigenvalues::Integer = 5,
                           k = nothing,
                           stencil::Integer = 9,
                           stencil_laplace::Integer = stencil,
                           stencil_nabla::Integer = stencil,
                           solver::Symbol = :arpack,
                           potential_unit::Unitful.Units = u"hartree",
                           coord_unit::Unitful.Units = u"bohr",
                           mass_unit::Unitful.Units = u"m_e")

    potential, system, output, files = setup_problem(V, coords;
        mass, periodic, n_eigenvalues, stencil, stencil_laplace, stencil_nabla,
        solver, potential_unit, coord_unit, mass_unit,
        reciprocal = k !== nothing)

    if k === nothing
        k_internal = Tuple(zeros(potential.dimension))
    else
        length(k) == potential.dimension ||
            throw(ArgumentError("k must have one component per dimension, got $(length(k)) for a $(potential.dimension)D problem"))
        # user k in coord_unit⁻¹ -> internal mass-weighted reciprocal atomic units
        k_bohr = ustrip.(uconvert.(u"bohr^-1", collect(float.(k)) .* coord_unit^-1))
        k_internal = Tuple(k_bohr ./ sqrt.(potential.mass))
    end

    energies, states = solve_at_k(potential, system, output, files, k_internal)

    kpoint = k === nothing ? nothing : collect(float.(k))
    return SchrodingerSolution(to_unit(energies, potential_unit), states, kpoint)
end

"""
    band_structure(V, coords; n_kpoints, kwargs...) -> BandStructure

Compute the band structure of a periodic system along the path through the
high-symmetry points of the Brillouin zone, without reading or writing any
files.

Takes the same potential/coordinate arguments and keyword arguments as
[`solve_schrodinger`](@ref) (except `k`), plus:

- `n_kpoints::Integer`: number of k-points per path segment (must be > 1).

At least one dimension must be `periodic`.
"""
function band_structure(V::AbstractArray{<:Real}, coords;
                        n_kpoints::Integer,
                        mass = 1.0,
                        periodic = true,
                        n_eigenvalues::Integer = 5,
                        stencil::Integer = 9,
                        stencil_laplace::Integer = stencil,
                        stencil_nabla::Integer = stencil,
                        solver::Symbol = :arpack,
                        potential_unit::Unitful.Units = u"hartree",
                        coord_unit::Unitful.Units = u"bohr",
                        mass_unit::Unitful.Units = u"m_e")

    n_kpoints > 1 || throw(ArgumentError("n_kpoints has to be larger than 1!"))

    potential, system, output, files = setup_problem(V, coords;
        mass, periodic, n_eigenvalues, stencil, stencil_laplace, stencil_nabla,
        solver, potential_unit, coord_unit, mass_unit,
        reciprocal = true)

    potential.n_kpoints     = n_kpoints
    potential.bandStructure = true
    kpoints = generate_kpoints(potential)

    n_k      = length(kpoints)
    energies = Matrix{Float64}(undef, n_k, n_eigenvalues)
    kpoints_user = Vector{Vector{Float64}}(undef, n_k)

    for (i, k) in enumerate(kpoints)
        e, _ = solve_at_k(potential, system, output, files, k)
        energies[i, :] .= to_unit(e, potential_unit)
        # internal mass-weighted k -> coord_unit⁻¹, as in the file pipeline
        k_bohr = collect(k) .* sqrt.(potential.mass)
        kpoints_user[i] = ustrip.(uconvert.(coord_unit^-1, k_bohr .* u"bohr^-1"))
    end

    kpath = zeros(n_k)
    for i in 2:n_k
        kpath[i] = kpath[i-1] + norm(kpoints_user[i] - kpoints_user[i-1])
    end

    return BandStructure(kpoints_user, kpath, energies)
end

####################################################################
#                                                                  #
# internals                                                        #
#                                                                  #
####################################################################

"""
Validate the array inputs and build fully initialized Potential, System,
Output and Files structs, mirroring exactly what readInputFile + checkInput +
readPotential + setupSystem produce for the file pipeline (see those files
for the unit conventions). Also builds the stencil operator matrices.
"""
function setup_problem(V::AbstractArray{<:Real}, coords;
                       mass, periodic, n_eigenvalues,
                       stencil, stencil_laplace, stencil_nabla,
                       solver, potential_unit, coord_unit, mass_unit,
                       reciprocal::Bool)

    dimension = ndims(V)
    1 <= dimension <= 3 ||
        throw(ArgumentError("only 1D, 2D and 3D potentials are supported, got a $(dimension)D array"))

    axes_ = coords isa AbstractVector{<:Real} ? (coords,) : Tuple(coords)
    length(axes_) == dimension ||
        throw(ArgumentError("expected $(dimension) coordinate axes for a $(dimension)D potential, got $(length(axes_))"))

    for (d, ax) in enumerate(axes_)
        length(ax) == size(V, d) ||
            throw(ArgumentError("coordinate axis $d has $(length(ax)) points but the potential has $(size(V, d)) along that dimension"))
        length(ax) >= 2 ||
            throw(ArgumentError("coordinate axis $d needs at least 2 points"))
        diffs = diff(collect(float.(ax)))
        all(>(0), diffs) ||
            throw(ArgumentError("coordinate axis $d must be strictly ascending"))
        isapprox(minimum(diffs), maximum(diffs); rtol = 1.0e-8) ||
            throw(ArgumentError("coordinate axis $d must be equally spaced"))
    end

    all(isfinite, V) || throw(ArgumentError("the potential contains non-finite values"))

    haskey(SOLVER_NAMES, solver) ||
        throw(ArgumentError("unknown solver :$solver - valid options are :arpack, :krylov and :lu"))
    solver === :arpack && n_eigenvalues + 5 >= length(V) &&
        throw(ArgumentError("the arpack solver needs n_eigenvalues + 5 < number of grid points ($(length(V)))"))

    stencil         in (3, 5, 7, 9, 11, 13) || throw(ArgumentError("stencil has to be 3, 5, 7, 9, 11 or 13"))
    stencil_laplace in (3, 5, 7, 9, 11, 13) || throw(ArgumentError("stencil-laplace has to be 3, 5, 7, 9, 11 or 13"))
    stencil_nabla   in (3, 5, 7, 9, 11)     || throw(ArgumentError("stencil-nabla has to be 3, 5, 7, 9 or 11"))

    n_eigenvalues >= 1 || throw(ArgumentError("n_eigenvalues has to be at least 1"))

    mass_vec = broadcast_per_dim(float.(mass), dimension, "mass")
    all(>(0), mass_vec) || throw(ArgumentError("masses must be positive"))
    periodic_vec = broadcast_per_dim(Bool.(periodic), dimension, "periodic")

    reciprocal && !any(periodic_vec) &&
        throw(ArgumentError("k-points require at least one periodic dimension"))

    potential = Potential()
    potential.dimension     = dimension
    potential.potentialUnit = potential_unit
    potential.coordsUnit    = coord_unit
    potential.massUnit      = mass_unit
    potential.mass          = ustrip.(uconvert.(u"m_e", mass_vec .* mass_unit))
    potential.periodic      = periodic_vec
    potential.n_datapoints  = collect(size(V))
    potential.reciprocal    = reciprocal
    potential.bandStructure = false
    potential.n_kpoints     = -1

    # flatten with the LAST dimension varying fastest (the file pipeline's
    # row-major convention, see buildStencilMatrices.jl) and convert units
    V_flat = dimension == 1 ? vec(V) : vec(permutedims(V, dimension:-1:1))
    potential.potential = ustrip.(uconvert.(u"hartree", float.(V_flat) .* potential_unit))

    axes_bohr = [ustrip.(uconvert.(u"bohr", collect(float.(ax)) .* coord_unit)) for ax in axes_]
    potential.coords = [expand_axis(axes_bohr, d) for d in 1:dimension]

    # mass-weighted grid spacing, readPotential.jl:103
    potential.intervall = [(axes_bohr[d][2] - axes_bohr[d][1]) * sqrt(potential.mass[d]) for d in 1:dimension]
    iv = potential.intervall
    all(isapprox.(iv, iv[1]; rtol = 1.0e-8)) ||
        throw(ArgumentError("the mass-weighted grid spacing must be equal in every dimension " *
                            "(step(coords[i]) * sqrt(mass[i]) is $(iv)); this is a restriction " *
                            "of the underlying discretization"))

    potential.shift = minimum(potential.potential)

    system = System()
    system.solver   = SOLVER_NAMES[solver]
    system.stencil  = stencil
    system.stencilΔ = stencil_laplace
    system.stencil∇ = stencil_nabla

    setupSystem(potential, system)

    buildΔ(system, potential)
    reciprocal && build∇(system, potential)

    output = Output()
    output.n_eigenvalues = n_eigenvalues

    files    = Files()
    files.to = TimerOutput()

    return potential, system, output, files
end

"""
Solve at one internal (mass-weighted) k-point, replicating the shift dance of
the main.jl k-loop, and return (energies in hartree with the potential shift
restored, states as normalized matrix columns).
"""
function solve_at_k(potential::Potential, system::System, output::Output, files::Files, k)
    potential.potential = potential.potential .- potential.shift
    try
        solve(potential, system, output, k, files)
    finally
        potential.potential = potential.potential .+ potential.shift
    end

    energies = output.eigenvalues .+ potential.shift
    states   = stack(output.eigenvectors)
    states   = system.reciprocal ? states : real.(states)

    # renormalize with the correct volume element: sum(abs2) * prod(Δx) == 1
    # with Δx per dimension in coord units (normalize.jl uses ψ² and a single
    # spacing, which is wrong for complex states and unequal plain spacings)
    spacings = potential.intervall ./ sqrt.(potential.mass)  # plain bohr
    volume_element = prod(ustrip(uconvert(potential.coordsUnit, dx * u"bohr")) for dx in spacings)
    for j in axes(states, 2)
        states[:, j] ./= sqrt(sum(abs2, view(states, :, j)) * volume_element)
    end

    return energies, states
end

"""
Generate the internal (mass-weighted) k-points for a band-structure run,
mirroring readPotential.jl:121-168.
"""
function generate_kpoints(potential::Potential)
    d = potential.dimension
    k_intervalls = [π / ((potential.coords[i][end] - potential.coords[i][1]) * sqrt(potential.mass[i]) +
                         potential.intervall[i]) / (potential.n_kpoints - 1) for i in 1:d]

    kpoints = d == 1 ? get_kpoints_1D(k_intervalls, potential.n_kpoints) :
              d == 2 ? get_kpoints_2D(k_intervalls, potential.n_kpoints) :
                       get_kpoints_3D(k_intervalls, potential.n_kpoints)

    # non-periodic dimensions only contribute k = 0, readPotential.jl:164-168
    for i in 1:d
        potential.periodic[i] || filter!(k -> k[i] == 0.0, kpoints)
    end

    isempty(kpoints) &&
        throw(ArgumentError("no k-points remain after removing non-periodic components"))

    return kpoints
end

"""
Per-point column of axis `d` for a grid flattened with the last dimension
varying fastest.
"""
function expand_axis(axes_bohr::Vector{Vector{Float64}}, d::Int)
    dims = Tuple(length.(axes_bohr))
    n = prod(dims)
    col = Vector{Float64}(undef, n)
    i = 0
    for ci in CartesianIndices(reverse(dims))
        idx = reverse(Tuple(ci))
        col[i += 1] = axes_bohr[d][idx[d]]
    end
    return col
end

function broadcast_per_dim(value, dimension::Int, name::String)
    vec = value isa AbstractVector ? collect(value) : fill(value, dimension)
    length(vec) == 1 && dimension > 1 && (vec = fill(vec[1], dimension))
    length(vec) == dimension ||
        throw(ArgumentError("$name must be a scalar or have one value per dimension, got $(length(vec)) values for a $(dimension)D problem"))
    return vec
end

to_unit(energies_hartree, unit) =
    ustrip.(uconvert.(unit, energies_hartree .* u"hartree"))

####################################################################
#                                                                  #
# Unitful convenience methods: quantities in, units inferred       #
#                                                                  #
####################################################################

"""
    solve_schrodinger(V::AbstractArray{<:Unitful.Energy}, coords; kwargs...)

Convenience method accepting a potential of `Unitful` energy quantities and
coordinate axes of length quantities; the units are inferred from the inputs.
"""
function solve_schrodinger(V::AbstractArray{<:Unitful.Energy}, coords; kwargs...)
    Vq = float.(ustrip.(u"hartree", V))
    solve_schrodinger(Vq, strip_length_axes(V, coords)...; potential_unit = u"hartree", coord_unit = u"bohr", kwargs...)
end

function band_structure(V::AbstractArray{<:Unitful.Energy}, coords; kwargs...)
    Vq = float.(ustrip.(u"hartree", V))
    band_structure(Vq, strip_length_axes(V, coords)...; potential_unit = u"hartree", coord_unit = u"bohr", kwargs...)
end

function strip_length_axes(V, coords)
    axes_ = coords isa AbstractVector ? (coords,) : Tuple(coords)
    stripped = map(axes_, ntuple(identity, length(axes_))) do ax, _
        eltype(ax) <: Unitful.Length ? ustrip.(u"bohr", ax) : ax
    end
    return (ndims(V) == 1 ? (stripped[1],) : (stripped,))
end
