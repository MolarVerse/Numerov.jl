function fileInfo(files, system)

    inputfile = files.inputFileName
    potfile   = files.PotentialFileName
    timefile  = files.timingsFileName

    stringbuffer = 
    "                                                                                        " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "       | input and output file information                                          |   " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "                                                                                        " * "\n" *
    "         input file    : \"$inputfile\"                                                 " * "\n" *
    "         potential file: \"$potfile\"                                                   " * "\n" *
    "                                                                                        " * "\n" *
    "         timings file  : \"$timefile\"                                                  " * "\n" *
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end