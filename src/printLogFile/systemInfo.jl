function systemInfo(files::Files, potential::Potential, system::System)
    
    dim           = potential.dimension
   
    coordsUnit    = potential.coordsUnit
    massUnit      = potential.massUnit
    potentialUnit = potential.potentialUnit

    stringbuffer = 
    "                                                                                        " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "       | System Information                                                         |   " * "\n" *
    "       ------------------------------------------------------------------------------   " * "\n" *
    "                                                                                        " * "\n" *
    "         dimension      = $(dim)D                                                       " * "\n" *
    "                                                                                        " * "\n" *
    "         coords unit    = $(coordsUnit)                                                 " * "\n" *
    "         mass unit      = $(massUnit)                                                   " * "\n" *
    "         potential unit = $(potentialUnit)                                              " * "\n" *
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end