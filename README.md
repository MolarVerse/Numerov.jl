# Numerov.jl

[![CI](https://github.com/MolarVerse/Numerov.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MolarVerse/Numerov.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/MolarVerse/Numerov.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MolarVerse/Numerov.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Numerov.jl solves the time-independent Schrödinger equation on 1D, 2D and 3D grid
potentials. The Hamiltonian is discretized with high-order finite-difference
(Numerov-type) stencils, assembled as sparse matrices and diagonalized with iterative
eigensolvers (Arpack, KrylovKit) or a dense LU-based solver. Periodic systems are
supported, including k-point sampling and band-structure calculations along paths
through the high-symmetry points of the Brillouin zone. The package was developed in
an academic computational-chemistry context, e.g. for vibrational eigenstates of
molecules on potential-energy surfaces from electronic-structure calculations.

## Requirements

- Julia >= 1.10

## Installation

The package is not yet registered in the Julia General registry (registration is
planned). Until then, install it directly from GitHub:

```julia-repl
pkg> add https://github.com/MolarVerse/Numerov.jl
```

## Quickstart

Clone the repository and instantiate its environment:

```sh
git clone https://github.com/MolarVerse/Numerov.jl.git
cd Numerov.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then run one of the ready-made examples:

```sh
cd examples/1DHarmonicOscillator
julia --project=../.. -e 'using Numerov; numerov("input.in")'
```

All output files are written to the current working directory. To keep the repository
clean, you can instead copy the example directory somewhere else, `cd` into the copy
and run the same command with `--project=<path to the cloned Numerov.jl>`.

The package exports a single function, `numerov(inputFileName::String)`, which reads
the input file, solves the Schrödinger equation and writes all result files (see
[Output files](#output-files)). Alternatively, a small launcher script is provided:

```sh
julia --project=<Numerov.jl dir> <Numerov.jl dir>/bin/numerov.jl input.in
```

Invalid input does not kill the Julia session: `numerov` throws an `ArgumentError`
that can be caught, so it is safe to call from the REPL or from your own scripts.

## Input file

The input file consists of `keyword = value` lines. Everything after a `#` is treated
as a comment and blank lines are ignored. Keywords are case-insensitive; values of
enumerated options (units, solver names, on/off switches) are also compared
case-insensitively, but file-name values are taken verbatim and are therefore
case-sensitive. Defining the same keyword twice is an error.

The only required keyword is `potential-file`. All others have defaults:

| Keyword | Meaning | Allowed values | Default |
| --- | --- | --- | --- |
| `potential-file` | Path to the grid-potential file | file path | **required** |
| `potential-unit` | Energy unit of the potential values | `hartree`, `ev`, `kj/mol`, `kcal/mol` | `hartree` |
| `coord-unit` | Length unit of the grid coordinates | `angstrom`, `nm`, `bohr` | `angstrom` |
| `mass-unit` | Unit of the reduced masses (`unit` and `g/mol` both mean the atomic mass unit, `me`/`m_e` is the electron mass) | `unit`, `g/mol`, `me`, `m_e` | `unit` |
| `stencil` | Number of finite-difference stencil points [^1] | `3`, `5`, `7`, `9`, `11`, `13` | `9` |
| `stencil-laplace` | Stencil size for the Laplacian only (overrides `stencil`) [^1] | `3`, `5`, `7`, `9`, `11`, `13` | value of `stencil` |
| `stencil-nabla` | Stencil size for the gradient only (overrides `stencil`; the gradient matrix is used for periodic/k-point runs) | `3`, `5`, `7`, `9`, `11` | value of `stencil` |
| `reduced-mass` | Reduced mass per dimension, comma- or space-separated; a single value is applied to all dimensions | real numbers | `1.0` |
| `periodic` | Periodic boundary conditions per dimension, comma- or space-separated; a single value is applied to all dimensions | `true`, `false` | `false` |
| `n-eigenvalues` | Number of eigenvalues/eigenvectors to compute | integer >= 1 | `5` |
| `k-points` | Number of k-points sampled per direction between the Gamma point and the Brillouin-zone boundary; if omitted, a single calculation at k = 0 is done | integer > 1 | not set |
| `datapoints` | Number of grid points per dimension, comma- or space-separated (e.g. `datapoints = 20, 30`) | integers | required for 2D/3D; in 1D taken from the potential file |
| `band-structure` | Compute the band structure along the path through the high-symmetry points of the Brillouin zone (requires `k-points`) | `on`, `true`, `off`, `false` | `off` |
| `solver` | Eigensolver backend | `arpack`, `krylov`, `lu` | `arpack` |
| `output-file` | Name of the log file | file path | `Numerov.out` |
| `timings-file` | Name of the timings file | file path | `timings.out` |
| `read-k-points` | Read the k-points from a file instead of generating them. **Accepted but not yet implemented** | `true`, `false` | `false` |
| `k-points-file` | k-point file for `read-k-points`. **Accepted but not yet implemented** -- the file is ignored with a warning | file path | not set |

[^1]: Not every stencil size is implemented for every dimensionality: the
    13-point Laplacian is not available in 2D or 3D, and the 3-point Laplacian
    is not available in 3D. Unsupported combinations throw an `ArgumentError`
    when the operator matrix is built.

An annotated example (based on `examples/1DKronigPenney/input.in`):

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

## Potential file format

The potential is given on a regular grid as whitespace-separated columns, one grid
point per line: first the coordinates, then the potential value. `#` starts a comment
and blank lines are ignored. The dimensionality of the problem is inferred from the
number of columns:

```text
# 1D:  x  V(x)
-10.0   50.0
 -9.9   49.005
 ...

# 2D:  x  y  V(x,y)
# 3D:  x  y  z  V(x,y,z)
```

Coordinates are interpreted in `coord-unit` and potential values in `potential-unit`;
internally everything is converted to atomic units. For 2D and 3D potentials the
number of grid points per dimension must be given explicitly with the `datapoints`
keyword.

## Output files

All output files are written to the current working directory:

| File | Content | Written |
| --- | --- | --- |
| `Numerov.out` | Log file: input parsing, system and sparse-matrix information (name configurable via `output-file`) | always |
| `eigenvalues.dat` | Eigenvalues in the chosen `potential-unit`, one line per k-point (k-point columns first for periodic runs) | always |
| `eigenvectors.dat`, `eigenvectors_shifted.dat` | Coordinates, potential and eigenvector amplitudes per grid point; the `shifted` variant offsets each eigenvector by its eigenvalue for plotting inside the potential | non-periodic runs |
| `eigenvectors_k_<k>.dat`, `eigenvectors_shifted_k_<k>.dat`, `imag_eigenvectors_k_<k>.dat`, `imag_eigenvectors_shifted_k_<k>.dat` | Real and imaginary parts of the (complex) eigenvectors, one set of files per k-point | periodic/k-point runs |
| `frequencies.dat` (or `frequencies_k_<k>.dat`) | Transition energies between the computed states as a lower triangular matrix, in cm^-1 | always (per k-point for periodic runs) |
| `bandstructure.dat` | Distance along the k-path and the corresponding eigenvalues | only with `band-structure = on` |
| `timings.out` | Timing breakdown of the calculation (name configurable via `timings-file`) | always |

Note that an existing `eigenvalues.dat` in the working directory is deleted at the
start of a run and replaced; other output files are overwritten as well.

## Examples

The `examples/` directory contains ready-to-run cases, each consisting of an
`input.in` and the corresponding potential data, including
`1DHarmonicOscillator`, `1DKronigPenney`, `2DKronigPenney`, `2DWater` (vibrational
states of water on a CCSD(T) surface), `3DKronigPenney` and more. The
`examples/generators/` directory provides `generate_cosine.jl` and
`generate_KronigPenney.jl` to build model potentials in 1D, 2D and 3D.

## Running tests

```julia-repl
pkg> test
```

This runs the unit tests and the 1D/2D integration tests. The two heavy 3D test cases
(3D harmonic oscillator and 3D Kronig-Penney) are skipped by default; set the
environment variable `NUMEROV_TEST_FULL=true` to include them.

## License

Numerov.jl is released under the [MIT License](LICENSE).
