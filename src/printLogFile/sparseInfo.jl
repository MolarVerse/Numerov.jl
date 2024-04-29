function sparseInfo(files::Files, system::System)
    
    n_non_zeros = length(system.Δ.nzval)
    n_zeros     = length(system.Δ) - length(system.Δ.nzval)
    sparsity    = length(system.Δ.nzval) / length(system.Δ) * 100

    stringbuffer = 
    "                                                                                        " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "       | Sparsity of Hamiltonian                                                    |   " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "                                                                                        " * "\n" *
    "         non zero entires = $n_non_zeros                                                " * "\n" *
    "         zero entries     = $n_zeros                                                    " * "\n" *
    "         sparsity         = $(@sprintf("%.3g%%", sparsity))                             " * "\n" *
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end