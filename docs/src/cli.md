# Command-line interface

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

The CLI writes the same output files as [`numerov`](@ref), into the current
working directory. Invalid input prints a single-line error message and exits
with a nonzero status instead of a stacktrace; when calling the library
function [`numerov`](@ref) from Julia instead, the same errors are thrown as
catchable `ArgumentError`s.

Without installing the shim, the same entry point is available as a script:

```sh
julia --project=<Numerov.jl dir> <Numerov.jl dir>/bin/numerov.jl input.in
```

## Entry point

```@docs
Numerov.CLI
Numerov.CLI.numerov
```
