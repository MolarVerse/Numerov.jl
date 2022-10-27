module Numerov

myexit() = run(`julia.exe`)

include("readInputFile.jl")
include("main.jl")

end # module Numerov
