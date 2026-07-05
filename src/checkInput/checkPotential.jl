function checkInput(potential::Potential)
    checkPotentialUnits(potential)
    checkCoordsUnits(potential)
    checkMassUnits(potential)
    checkMass(potential)
    checkKPoints(potential)
    checkNDatapoints(potential)
    checkBandStructure(potential)
    checkPeriodicity(potential)
    checkReadKPoints(potential)
end

function checkPotentialUnits(potential::Potential)
    potentialUnit = lowercase(inputDictionary["potential-unit"])
    isempty(potentialUnit)      && (potential.potentialUnit = u"hartree"    ; return) #write to log file about default setting
    potentialUnit == "hartree"  && (potential.potentialUnit = u"hartree"    ; return)
    potentialUnit == "ev"       && (potential.potentialUnit = u"eV"         ; return)
    potentialUnit == "kj/mol"   && (potential.potentialUnit = u"kJpermol"   ; return)
    potentialUnit == "kcal/mol" && (potential.potentialUnit = u"kcalpermol" ; return)

    throw(ArgumentError("\nThe given potential-unit $(inputDictionary["potential-unit"]) was not recognised!\n" *
                        "Valid options are:                                                                 \n" *
                        "    - hartree                                                                      \n" *
                        "    - ev                                                                           \n" *
                        "    - kj/mol                                                                       \n" *
                        "    - kcal/mol                                                                     \n"))
end

function checkCoordsUnits(potential::Potential)
    coordsUnit = lowercase(inputDictionary["coord-unit"])
    isempty(coordsUnit)      && (potential.coordsUnit = u"angstrom"; return) #write to log file about default setting
    coordsUnit == "angstrom" && (potential.coordsUnit = u"angstrom"; return)
    coordsUnit == "nm"       && (potential.coordsUnit = u"nm"      ; return)
    coordsUnit == "bohr"     && (potential.coordsUnit = u"bohr"    ; return)

    throw(ArgumentError("\nThe given coord-unit $(inputDictionary["coord-unit"]) was not recognised!\n" *
                        "Valid options are:                                                         \n" *
                        "    - angstrom                                                             \n" *
                        "    - nm                                                                   \n" *
                        "    - bohr                                                                 \n"))
end

function checkMassUnits(potential::Potential)
    massUnit = lowercase(inputDictionary["mass-unit"])
    isempty(massUnit)   && (potential.massUnit = u"u" ; return) #write to log file about default setting
    massUnit == "unit"  && (potential.massUnit = u"u" ; return)
    massUnit == "g/mol" && (potential.massUnit = u"u" ; return)
    massUnit == "me"    && (potential.massUnit = u"me"; return)

    throw(ArgumentError("\nThe given mass-unit $(inputDictionary["mass-unit"]) was not recognised!\n" *
                        "Valid options are:                                                       \n" *
                        "    - unit                                                               \n" *
                        "    - g/mol                                                              \n" *
                        "    - me            #stands for electron mass                            \n"))
end

function checkMass(potential::Potential)
    isempty(inputDictionary["reduced-mass"]) && (potential.mass = [1.0]; return) #write to log file about default setting
    potential.mass = parse.(Float64, split(join(split(inputDictionary["reduced-mass"], ","), " ")))
    length(potential.mass) == 0 && throw(ArgumentError("reduced-mass not correctly defined!"))
end

function checkKPoints(potential::Potential)
    isempty(inputDictionary["k-points"]) && (potential.n_kpoints = -1; return) #write to log file about default setting
    potential.n_kpoints = parse(Int64, inputDictionary["k-points"])
    potential.n_kpoints < 2 && throw(ArgumentError("number of k-points has to be > 1"))
end

function checkNDatapoints(potential::Potential)
    isempty(inputDictionary["datapoints"]) && (potential.n_datapoints = Vector(); return) #write to log file about default setting
    potential.n_datapoints = parse.(Int64, split(join(split(inputDictionary["datapoints"], ","), " ")))
    length(potential.n_datapoints) == 0 && throw(ArgumentError("n_datapoints not correctly defined!"))
end

function checkBandStructure(potential::Potential)
    bandStructure = lowercase(inputDictionary["band-structure"])
    isempty(bandStructure)            && (potential.bandStructure = false; return) #write to log file about default setting
    bandStructure ∈ ["on" , "true" ]  && (potential.bandStructure = true ; return)
    bandStructure ∈ ["off", "false"]  && (potential.bandStructure = false; return)

    throw(ArgumentError("\nThe given band-structure option $(inputDictionary["band-structure"]) was not recognised!\n" *
                        "Valid options are: on, true, off, false"))
end

function checkPeriodicity(potential::Potential)
    isempty(inputDictionary["periodic"]) && (potential.periodic = [false]; return) #write to log file about default setting
    potential.periodic = parse.(Bool, split(join(split(lowercase(inputDictionary["periodic"]), ","), " ")))
end

function checkReadKPoints(potential::Potential)
    isempty(inputDictionary["read-k-points"]) && (potential.read_kpoints = false; return)
    potential.read_kpoints = parse(Bool, lowercase(inputDictionary["read-k-points"]))
    if potential.read_kpoints
        if !isempty(inputDictionary["k-points"])
            throw(ArgumentError("reading from k-points file is only allowed if no k-points are given!"))
        elseif isempty(inputDictionary["k-points-file"])
            throw(ArgumentError("to choose option to read from k-points file you have to specify the file in the input file via the keyword \"k-points-file\"!"))
        end
    end
end
