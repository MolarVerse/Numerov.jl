# Library usage

Besides the input-file pipeline around [`numerov`](@ref), Numerov.jl can be
used as a plain Julia library: [`solve_schrodinger`](@ref) and
[`band_structure`](@ref) take the potential as an array and the grid axes as
ranges or vectors, and return the results as ordinary Julia values. No files
are read or written and no global state is touched, so the functions are safe
to call repeatedly, e.g. inside a parameter scan.

## A first example: the 1D harmonic oscillator

For the harmonic potential ``V(x) = \tfrac{1}{2} x^2`` in atomic units the
exact eigenvalues are ``E_n = n + \tfrac{1}{2}``:

```julia
using Numerov

x = range(-5.0, 5.0; length = 201)     # grid in bohr (the default coord unit)
V = 0.5 .* x .^ 2                      # potential in hartree (the default)

result = solve_schrodinger(V, x; n_eigenvalues = 3)

result.energies   # ≈ [0.5, 1.5, 2.5] hartree
result.states     # 201×3 matrix, one normalized eigenvector per column
```

The returned [`SchrodingerSolution`](@ref) holds the eigenvalues (`energies`),
the eigenvectors as matrix columns (`states`) and, for periodic calculations,
the k-point the equation was solved at (`kpoint`).

## Grids and dimensionality

The potential is a 1-, 2- or 3-dimensional array and `coords` provides one
axis per array dimension: `V[i, j, k]` is the value at
`(coords[1][i], coords[2][j], coords[3][k])`. Each axis must be strictly
ascending and equally spaced. In 1D the axis may be passed directly; for 2D
and 3D pass a tuple (or vector) of axes.

A 2D isotropic harmonic oscillator shows the expected degeneracies
``E_{n_x n_y} = n_x + n_y + 1``:

```julia
x = range(-5.0, 5.0; length = 101)
y = range(-5.0, 5.0; length = 101)
V = [0.5 * (xi^2 + yj^2) for xi in x, yj in y]

result = solve_schrodinger(V, (x, y); n_eigenvalues = 6)

result.energies   # ≈ [1.0, 2.0, 2.0, 3.0, 3.0, 3.0]
```

## Periodic systems and band structures

Setting `periodic = true` (per dimension) imposes periodic boundary
conditions. On a periodic grid the last point must not duplicate the first
one — the point after the last grid point *is* the first one again, one step
further:

```julia
n = 64
x = range(-π, π; length = n + 1)[1:n]   # duplicated endpoint removed
V = cos.(x) .+ 1
```

A single k-point is solved with the `k` keyword (components in inverse
coordinate units); the states are then complex Bloch states:

```julia
sol = solve_schrodinger(V, x; periodic = true, k = [0.25], n_eigenvalues = 3)
sol.kpoint            # [0.25]
eltype(sol.states)    # ComplexF64
```

[`band_structure`](@ref) samples the path through the high-symmetry points of
the Brillouin zone with `n_kpoints` points per path segment and returns a
[`BandStructure`](@ref) whose `energies` matrix has one row per k-point and
one column per band. `kpath` is the cumulative distance along the path — the
natural x axis of a band-structure plot:

```julia
bands = band_structure(V, x; n_kpoints = 30, n_eigenvalues = 4)

size(bands.energies)   # (30, 4)

using Plots
plot(bands.kpath, bands.energies;
     xlabel = "k path (bohr⁻¹)", ylabel = "E (hartree)", legend = false)
```

## Units

The arrays are plain numbers; their units are declared with keyword
arguments. The defaults are atomic units — `potential_unit = u"hartree"`,
`coord_unit = u"bohr"`, `mass_unit = u"m_e"` (electron masses) — and the
energies are returned in `potential_unit`:

```julia
using Numerov
using Unitful

x = range(-2.0, 2.0; length = 201)     # angstrom
V = 100.0 .* x .^ 2                    # eV

result = solve_schrodinger(V, x;
                           potential_unit = u"eV",
                           coord_unit    = u"angstrom",
                           mass_unit     = u"u",       # atomic mass unit
                           mass          = 1.0)

result.energies   # in eV
```

`u"eV"`, `u"angstrom"` and `u"u"` ship with
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl); the atomic units
`u"bohr"`, `u"hartree"` and `u"me_au"` additionally require
`using UnitfulAtomic` in your code.

Alternatively, pass Unitful quantities directly — a potential of energy
quantities and axes of length quantities. The units are inferred from the
inputs, everything is converted to atomic units, and the energies come back
in hartree:

```julia
using Numerov
using Unitful

x = collect(range(-2.0, 2.0; length = 201)) .* u"angstrom"
V = 100.0 .* range(-2.0, 2.0; length = 201) .^ 2 .* u"eV"

result = solve_schrodinger(V, x; mass = 1.0, mass_unit = u"u")

result.energies                              # in hartree
uconvert.(u"eV", result.energies .* u"hartree")
```

## Normalization of the states

Each column of `states` is normalized so that

```julia
sum(abs2, state) * prod(spacings) == 1
```

with one grid spacing per dimension, expressed in `coord_unit`. This is the
discrete form of ``\int |\psi|^2 \, dV = 1`` on the grid. States are real for
non-periodic calculations and complex for periodic ones.

## The equal mass-weighted spacing restriction

The underlying discretization works on a mass-weighted grid and requires the
same effective spacing in every dimension:
`step(coords[i]) * sqrt(mass[i])` must be equal for all `i`, otherwise an
`ArgumentError` is thrown. Unequal spacings are therefore only possible when
the masses compensate them:

```julia
x = range(-5.0, 5.0; length = 101)     # step 0.1
y = range(-2.5, 2.5; length = 101)     # step 0.05
V = [0.5 * (xi^2 + 4 * yj^2) for xi in x, yj in y]

# 0.1 * sqrt(1.0) == 0.05 * sqrt(4.0) — accepted
result = solve_schrodinger(V, (x, y); mass = [1.0, 4.0], n_eigenvalues = 3)
```

For equal masses, simply use the same grid spacing in every dimension.

## When to still use the input-file pipeline

[`numerov`](@ref) and the [command-line interface](cli.md) read an
[input file](input.md) and write the full set of result files —
`eigenvalues.dat`, plot-ready (shifted) eigenvector files, transition
frequencies in cm⁻¹ and `bandstructure.dat`. Prefer that pipeline when
running from the shell, when the results should end up in files anyway, or
when a calculation is defined by a small self-documenting input file that can
be archived and re-run reproducibly. Prefer the programmatic API when the
potential is generated in Julia, when the results feed directly into further
Julia code, or when many related problems are solved in a loop.
