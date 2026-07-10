using Documenter
using Numerov

makedocs(
    sitename = "Numerov.jl",
    modules  = [Numerov],
    pages    = [
        "Home"                   => "index.md",
        "Input file reference"   => "input.md",
        "Command-line interface" => "cli.md",
        "API reference"          => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://molarverse.github.io/Numerov.jl",
    ),
    checkdocs = :exports,
    warnonly  = [:missing_docs],
)

deploydocs(
    repo      = "github.com/MolarVerse/Numerov.jl.git",
    devbranch = "main",
)
