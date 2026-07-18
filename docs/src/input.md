# Input file reference

The input file consists of `keyword = value` lines. `#` starts a comment and
blank lines are ignored. Keywords and enumerated values (units, solver names,
on/off switches) are case-insensitive. File-name values are the exception: they
are taken verbatim and are therefore case-sensitive. Defining the same keyword
twice is an error.

A typical input file (from `examples/1DKronigPenney/input.in`, annotated):

```text
# 1D Kronig-Penney model, periodic, with k-point sampling
potential-file = potential.dat     # grid potential (required)
stencil        = 9                 # 9-point finite-difference stencil
reduced-mass   = 23.06000881918738 # in units of mass-unit
mass-unit      = me                # electron masses
coord-unit     = angstrom          # unit of the coordinates in potential.dat
potential-unit = ev                # unit of the potential values in potential.dat
n-eigenvalues  = 10                # number of states to compute
periodic       = true              # periodic boundary conditions
k-points       = 10                # k-points from Gamma to the zone boundary
```

## Keywords

The only required keyword is `potential-file`. All others have defaults:

| Keyword | Meaning | Allowed values | Default |
| --- | --- | --- | --- |
| `potential-file` | Path to the grid-potential file | file path | **required** |
| `potential-unit` | Energy unit of the potential values | `hartree`, `ev`, `kj/mol`, `kcal/mol` | `hartree` |
| `coord-unit` | Length unit of the grid coordinates | `angstrom`, `nm`, `bohr` | `angstrom` |
| `mass-unit` | Unit of the reduced masses (`unit` and `g/mol` both mean the atomic mass unit; `me` and `m_e` both mean the electron mass) | `unit`, `g/mol`, `me`, `m_e` | `unit` |
| `stencil` | Number of finite-difference stencil points (see note below) | `3`, `5`, `7`, `9`, `11`, `13` | `9` |
| `stencil-laplace` | Stencil size for the Laplacian only (overrides `stencil`) | `3`, `5`, `7`, `9`, `11`, `13` | value of `stencil` |
| `stencil-nabla` | Stencil size for the gradient only (overrides `stencil`; the gradient matrix is used for periodic/k-point runs) | `3`, `5`, `7`, `9`, `11` | value of `stencil` |
| `reduced-mass` | Reduced mass per dimension, comma- or space-separated; a single value is applied to all dimensions | real numbers | `1.0` |
| `periodic` | Periodic boundary conditions per dimension, comma- or space-separated; a single value is applied to all dimensions | `true`, `false` | `false` |
| `n-eigenvalues` | Number of eigenvalues/eigenvectors to compute | integer >= 1 | `5` |
| `k-points` | Number of k-points sampled per direction between the Gamma point and the Brillouin-zone boundary; if omitted, a single calculation at k = 0 is performed | integer > 1 | not set |
| `datapoints` | Number of grid points per dimension, comma- or space-separated (e.g. `datapoints = 20, 30`) | integers | required for 2D/3D; in 1D taken from the potential file |
| `band-structure` | Compute the band structure along the path through the high-symmetry points of the Brillouin zone (requires `k-points`) | `on`, `true`, `off`, `false` | `off` |
| `solver` | Eigensolver backend (`lobpcg` is recommended for large non-periodic 3D problems) | `arpack`, `krylov`, `lobpcg`, `lu` | `arpack` |
| `output-file` | Name of the log file | file path | `Numerov.out` |
| `timings-file` | Name of the timings file | file path | `timings.out` |
| `read-k-points` | Read the k-points from a file instead of generating them (accepted but not yet implemented) | `true`, `false` | `false` |
| `k-points-file` | k-point file for `read-k-points` (accepted but not yet implemented — the file is ignored with a warning) | file path | not set |

!!! note "Stencil sizes and dimensionality"
    Not every stencil size is implemented for every dimensionality: the
    13-point Laplacian is not available in 2D or 3D, and the 3-point Laplacian
    is not available in 3D. Unsupported combinations throw an `ArgumentError`
    when the operator matrix is built.

## Potential file format

The potential is given on a regular grid as whitespace-separated columns, one
grid point per line: first the coordinates, then the potential value. `#`
starts a comment and blank lines are ignored. The dimensionality of the problem
is inferred from the number of columns:

```text
# 1D:  x  V(x)
-10.0   50.0
 -9.9   49.005
 ...

# 2D:  x  y  V(x,y)
# 3D:  x  y  z  V(x,y,z)
```

Coordinates are interpreted in `coord-unit` and potential values in
`potential-unit`; internally everything is converted to atomic units. For 2D
and 3D potentials the number of grid points per dimension must be given
explicitly with the `datapoints` keyword.
