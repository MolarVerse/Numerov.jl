# Numerov.jl

[![version](https://juliahub.com/docs/General/Numerov/stable/version.svg)](https://juliahub.com/ui/Packages/General/Numerov)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21389968.svg)](https://doi.org/10.5281/zenodo.21389968)
[![docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://molarverse.github.io/Numerov.jl/stable/)
[![docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://molarverse.github.io/Numerov.jl/dev/)
[![CI](https://github.com/MolarVerse/Numerov.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/MolarVerse/Numerov.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/MolarVerse/Numerov.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MolarVerse/Numerov.jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Numerov.jl solves the time-independent Schrödinger equation on 1D, 2D and 3D grid
potentials — for example, vibrational eigenstates of molecules on potential-energy
surfaces from electronic-structure calculations, or periodic model systems with
k-point sampling and band structures along paths through the high-symmetry points
of the Brillouin zone.

Under the hood, the Hamiltonian is discretized with high-order finite-difference
(Numerov-type) stencils, assembled as a sparse matrix and diagonalized with
iterative eigensolvers (Arpack, KrylovKit) or a dense LU-based solver.

![Harmonic-oscillator eigenstates and Kronig-Penney band structure computed with Numerov.jl](docs/assets/showcase.png)

*Left: the five lowest eigenstates of a 1D harmonic oscillator, computed from
`examples/1DHarmonicOscillator`. Right: the band structure of a 1D Kronig-Penney
model from `examples/1DKronigPenney`.*

## Installation

Requires Julia 1.10 or later. Numerov.jl is registered in the Julia General
registry:

```julia-repl
pkg> add Numerov
```

Installing this way is all you need to use the package. The ready-to-run
examples ship with the repository, so to try one of those, clone the repository
instead (or in addition) — see the Quickstart below.

## Quickstart

Numerov.jl has a single entry point: the exported function
`numerov(inputFileName::String)` reads an input file, solves the Schrödinger
equation and writes all result files (see [Output files](#output-files)) to the
current working directory. On invalid input it throws a catchable
`ArgumentError` rather than terminating Julia, so it is safe to call from the
REPL or from your own scripts.

To run a bundled example, clone the repository and instantiate its environment:

```sh
git clone https://github.com/MolarVerse/Numerov.jl.git
cd Numerov.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then solve the 1D harmonic oscillator:

```sh
cd examples/1DHarmonicOscillator
julia --project=../.. -e 'using Numerov; numerov("input.in")'
```

This writes `eigenvalues.dat` with the first ten eigenvalues (0.5, 1.5, 2.5, …
E<sub>h</sub> — the exact harmonic-oscillator spectrum), `eigenvectors_shifted.dat`
for plotting the states inside the well (as in the figure above), and the other
files listed under [Output files](#output-files).

> **Tip:** output files land in the current working directory. To keep the
> repository clean, copy an example directory somewhere else and run the same
> command with `--project=<Numerov.jl dir>`.

A small launcher script is also provided:

```sh
julia --project=<Numerov.jl dir> <Numerov.jl dir>/bin/numerov.jl input.in
```

## Library usage

The Schrödinger equation can also be solved directly from Julia arrays,
without input or output files:

```julia
using Numerov

x = range(-5.0, 5.0; length = 201)   # grid in bohr
result = solve_schrodinger(0.5 .* x .^ 2, x; n_eigenvalues = 3)

result.energies   # ≈ [0.5, 1.5, 2.5] hartree
result.states     # normalized eigenvectors, one per column
```

`band_structure` computes band structures of periodic potentials the same
way. See the [library usage documentation](https://molarverse.github.io/Numerov.jl/dev/library/)
for 2D/3D problems, periodic systems, units and the full keyword list.

## Command-line interface

Numerov.jl ships a [Comonicon](https://github.com/comonicon/Comonicon.jl)-based
command-line interface. Install the `numerov` command once:

```julia-repl
julia> using Numerov
julia> Numerov.CLI.comonicon_install()
```

This places a `numerov` executable in `~/.julia/bin` (add that directory to
your `PATH`). Then:

```sh
numerov input.in      # run a calculation
numerov --help        # usage and argument description
numerov --version     # package version
```

Invalid input prints a single-line error message and exits with a nonzero
status instead of a stacktrace.

## Input file

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

The only required keyword is `potential-file`. All others have defaults:

| Keyword | Meaning | Allowed values | Default |
| --- | --- | --- | --- |
| `potential-file` | Path to the grid-potential file | file path | **required** |
| `potential-unit` | Energy unit of the potential values | `hartree`, `ev`, `kj/mol`, `kcal/mol` | `hartree` |
| `coord-unit` | Length unit of the grid coordinates | `angstrom`, `nm`, `bohr` | `angstrom` |
| `mass-unit` | Unit of the reduced masses [^1] | `unit`, `g/mol`, `me`, `m_e` | `unit` |
| `stencil` | Number of finite-difference stencil points [^2] | `3`, `5`, `7`, `9`, `11`, `13` | `9` |
| `stencil-laplace` | Stencil size for the Laplacian only (overrides `stencil`) [^2] | `3`, `5`, `7`, `9`, `11`, `13` | value of `stencil` |
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
| `read-k-points` | Read the k-points from a file instead of generating them [^3] | `true`, `false` | `false` |
| `k-points-file` | k-point file for `read-k-points` [^3] | file path | not set |

[^1]: `unit` and `g/mol` both mean the atomic mass unit; `me` and `m_e` both
    mean the electron mass.

[^2]: Not every stencil size is implemented for every dimensionality: the
    13-point Laplacian is not available in 2D or 3D, and the 3-point Laplacian
    is not available in 3D. Unsupported combinations throw an `ArgumentError`
    when the operator matrix is built.

[^3]: Accepted but not yet implemented — the file is ignored with a warning.

## Potential file format

The potential is given on a regular grid as whitespace-separated columns, one grid
point per line: first the coordinates, then the potential value. `#` starts a
comment and blank lines are ignored. The dimensionality of the problem is inferred
from the number of columns:

```text
# 1D:  x  V(x)
-10.0   50.0
 -9.9   49.005
 ...

# 2D:  x  y  V(x,y)
# 3D:  x  y  z  V(x,y,z)
```

Coordinates are interpreted in `coord-unit` and potential values in
`potential-unit`; internally everything is converted to atomic units. For 2D and
3D potentials the number of grid points per dimension must be given explicitly
with the `datapoints` keyword.

## Output files

All output files are written to the current working directory; existing files
with the same names are overwritten on each run (`eigenvalues.dat` is deleted as
soon as a run starts, since it is appended to per k-point).

| File | Content | Written |
| --- | --- | --- |
| `Numerov.out` | Log file: input parsing, system and sparse-matrix information (name configurable via `output-file`) | always |
| `eigenvalues.dat` | Eigenvalues in the chosen `potential-unit`, one line per k-point (k-point columns first for periodic runs) | always |
| `eigenvectors.dat`, `eigenvectors_shifted.dat` | Coordinates, potential and eigenvector amplitudes per grid point; the `shifted` variant offsets each eigenvector by its eigenvalue for plotting inside the potential | non-periodic runs |
| `eigenvectors_k_<k>.dat`, `eigenvectors_shifted_k_<k>.dat`, `imag_eigenvectors_k_<k>.dat`, `imag_eigenvectors_shifted_k_<k>.dat` | Real and imaginary parts of the (complex) eigenvectors, one set of files per k-point | periodic/k-point runs |
| `frequencies.dat` (or `frequencies_k_<k>.dat`) | Transition energies between the computed states as a lower triangular matrix, in `cm⁻¹` | always (per k-point for periodic runs) |
| `bandstructure.dat` | Distance along the k-path and the corresponding eigenvalues | only with `band-structure = on` |
| `timings.out` | Timing breakdown of the calculation (name configurable via `timings-file`) | always |

## Examples

The `examples/` directory contains ready-to-run cases, each with an `input.in`
and its potential data. Highlights: `1DHarmonicOscillator`, `1DKronigPenney`,
`2DKronigPenney`, `2DWater` (vibrational states of water on a CCSD(T) surface)
and `3DKronigPenney`. The `examples/generators/` directory provides
`generate_cosine.jl` and `generate_KronigPenney.jl` to build model potentials in
1D, 2D and 3D.

## Running tests

From the cloned repository, start `julia --project=.` and run:

```julia-repl
pkg> test
```

This runs the unit tests and the 1D/2D integration tests. The two heavy 3D test
cases (3D harmonic oscillator and 3D Kronig-Penney) are skipped by default; set
the environment variable `NUMEROV_TEST_FULL=true` to include them.

## Citing

If you use Numerov.jl in your research, please cite:

> J. Gamper, F. Kluibenschedl, A. K. H. Weiss, T. S. Hofer,
> *Accessing Position Space Wave Functions in Band Structure Calculations of
> Periodic Systems — A Generalized, Adapted Numerov Implementation for One-,
> Two-, and Three-Dimensional Quantum Problems*,
> J. Phys. Chem. Lett. **2023**, 14, 33, 7395–7403.
> [doi:10.1021/acs.jpclett.3c01707](https://doi.org/10.1021/acs.jpclett.3c01707)

Citation metadata is also provided in [`CITATION.cff`](CITATION.cff) (GitHub's
"Cite this repository" button generates BibTeX from it).

## Contributing

Bug reports and pull requests are welcome on the
[issue tracker](https://github.com/MolarVerse/Numerov.jl/issues).

## License

Numerov.jl is released under the [MIT License](LICENSE).
