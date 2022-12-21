function fileInfo(files::Files, potential::Potential, system::System)

    inputfile = files.inputFileName
    potfile   = files.PotentialFileName

    outputfile = files.logFileName 
    timefile   = files.timingsFileName

    eigenvaluefile   = files.eigenvalueFileName
    if system.reciprocal
        eigenvectorfiles              = "eigenvectors_k_*.dat" 
        imag_eigenvectorfiles         = "imag_eigenvectors_k_*.dat" 
        shifted_eigenvectorfiles      = "eigenvectors_shifted_k_*.dat" 
        shifted_imag_eigenvectorfiles = "imag_eigenvectors_shifted_k_*.dat"
        frequencyfiles                = "frequencies_k_*.dat"
    else
        eigenvectorfiles              = "eigenvectors.dat" 
        imag_eigenvectorfiles         = "-"
        shifted_eigenvectorfiles      = "eigenvectors_shifted.dat" 
        shifted_imag_eigenvectorfiles = "-" 
        frequencyfiles                = "frequencies.dat"
    end

    if potential.bandStructure
        bandstructurefile = files.bandStructureFileName
    else
        bandstructurefile = "-"
    end

    stringbuffer = 
    "                                                                                        " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "       | input and output file information                                          |   " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "                                                                                        " * "\n" *
    "         #########################  INPUT FILES  ##################################     " * "\n" *
    "                                                                                        " * "\n" *
    "         input file                       : \"$inputfile\"                              " * "\n" *
    "         potential file                   : \"$potfile\"                                " * "\n" *
    "                                                                                        " * "\n" *
    "         #########################  OUTPUT FILES  #################################     " * "\n" *
    "                                                                                        " * "\n" *
    "         output file                      : \"$outputfile\"                             " * "\n" *
    "         timings file                     : \"$timefile\"                               " * "\n" *
    "                                                                                        " * "\n" *
    "         eigenvalue file                  : \"$eigenvaluefile\"                         " * "\n" *
    "                                                                                        " * "\n" *
    "         real eigenvector file(s)         : \"$eigenvectorfiles\"                       " * "\n" *
    "         imag eigenvector file(s)         : \"$imag_eigenvectorfiles\"                  " * "\n" *
    "                                                                                        " * "\n" *
    "         shifted real eigenvector file(s) : \"$shifted_eigenvectorfiles\"               " * "\n" *
    "         shifted imag eigenvector file(s) : \"$shifted_imag_eigenvectorfiles\"          " * "\n" *
    "                                                                                        " * "\n" *
    "         frequencies file(s)              : \"$frequencyfiles\"                         " * "\n" *
    "                                                                                        " * "\n" *
    "         bandstructure file               : \"$bandstructurefile\"                      " * "\n" *
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end