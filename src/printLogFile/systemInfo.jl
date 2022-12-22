function systemInfo(files::Files, potential::Potential, system::System)
    
    dim           = potential.dimension
   
    coordsUnit    = potential.coordsUnit
    massUnit      = potential.massUnit
    potentialUnit = potential.potentialUnit

    a = map(x -> x ? "x" : "-", potential.periodic)

    while length(a) < 3
        push!(a, " ")
    end

    b = map(x -> x && potential.reciprocal ? "x" : "-", potential.periodic) #TODO: make reciprocal into all periodic directions

    while length(b) < 3
        push!(b, " ")
    end

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
    "                                                                                        " * "\n" *
    "                    |   x   |   y   |   z   |                                           " * "\n" *
    "                    -------------------------                                           " * "\n" *
    "         periodic   |   $(a[1])   |   $(a[2])   |   $(a[3])   |                         " * "\n" *
    "         reciprocal |   $(b[1])   |   $(b[2])   |   $(b[3])   |                         " * "\n" *
    "         n datapoints |   $(b[1])   |   $(b[2])   |   $(b[3])   |                         " * "\n" * #TODO: make this table dynamic!!!
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end