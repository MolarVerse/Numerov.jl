function systemInfo(files::Files, potential::Potential, system::System)
    
    dim           = potential.dimension
   
    coordsUnit    = potential.coordsUnit
    massUnit      = potential.massUnit
    potentialUnit = potential.potentialUnit

    xyz    = repeat([""], 5)
    xyz[2] = "x"
    xyz[3] = "y"
    xyz[4] = "z"

    max_length_string = calc_max_length_string(vcat(potential.n_datapoints, potential.n_kpoints))

    max_length_string = max_length_string %2 == 0 ? max_length_string + 1 : max_length_string

    header     = ["x","y","z"]
    delimiter  = "-"

    a = map(x -> x ? "x" : "-", potential.periodic)

    while length(a) < 3
        push!(a, " ")
    end

    b = map(x -> x && potential.reciprocal ? "x" : "-", potential.periodic) #TODO: make reciprocal into all periodic directions
    n = string.(potential.n_datapoints)
    k = [potential.periodic[i] ? string(potential.n_kpoints) : "0" for i in 1:dim]

    while length(b) < 3
        push!(b, " ")
    end
    while length(n) < 3
        push!(n, " ")
    end
    while length(k) < 3
        push!(k, " ")
    end

    for i in 1:max_length_string÷2
        header = map(x -> " " * x * " ", header)
        a = map(x -> " " * x * " ", a)
        b = map(x -> " " * x * " ", b)
        delimiter = "-" * delimiter * "-"
    end

    for i in eachindex(n)
        while length(n[i]) < max_length_string
            n[i] = " " * n[i]
        end
    end

    for i in eachindex(k)
        while length(k[i]) < max_length_string
            k[i] = " " * k[i]
        end
    end

    bandStructure = potential.bandStructure ? "yes" : "no"

    μ = Vector()
    for μ_i in potential.mass
        push!(μ, string(μ_i*uconvert(potential.massUnit, 1.0*potential.internalElemMass)))
    end

    while length(μ) < 3
        push!(μ, "-")
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
    "                      |   $(header[1])   |   $(header[2])   |   $(header[3])   |        " * "\n" *
    "                      ----$(delimiter)---|---$(delimiter)---|---$(delimiter)---|        " * "\n" *
    "         periodic     |   $(a[1])   |   $(a[2])   |   $(a[3])   |                       " * "\n" *
    "         reciprocal   |   $(b[1])   |   $(b[2])   |   $(b[3])   |                       " * "\n" *
    "         n-datapoints |   $(n[1])   |   $(n[2])   |   $(n[3])   |                       " * "\n" *
    "         n-kpoints    |   $(k[1])   |   $(k[2])   |   $(k[3])   |                       " * "\n" *
    "                                                                                        " * "\n" *
    "         band structure calculation requested: $(bandStructure)                         " * "\n" *
    "                                                                                        " * "\n" *
    "         μx = $(μ[1])                                                                   " * "\n" *
    "         μy = $(μ[2])                                                                   " * "\n" *
    "         μz = $(μ[3])                                                                   " * "\n" *
    "                                                                                        " * "\n"

    print(files.logFile, stringbuffer)
end

function determine_header()
    
end

calc_max_length_string(a) = maximum(length.(string.(a)))