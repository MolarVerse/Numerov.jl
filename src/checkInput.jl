function checkInput(potential::Potential)
    checkPotentialFile(potential)
    checkPotentialUnits(potential)
    checkCoordsUnits(potential)
    checkMassUnits(potential)
    checkMass(potential)
    checkShiftPotential(potential)
    checkKPoints(potential)
    checkNDatapoints(potential)
    checkBandStructure(potential)
end

function checkInput(system::System)
    checkStencil(system)
    checkPeriodicity(system)
end

function checkInput(output::Output)
    checkNEigenvalues(output)
end

function checkPotentialFile(potential::Potential)
    isempty(inputDictionary["potential-file"]) && (@error "No potential input file name given"; exit())
    potential.file = inputDictionary["potential-file"]
end

function checkPotentialUnits(potential::Potential)
    isempty(inputDictionary["potential-unit"])      && (potential.potentialUnit = u"hartree"    ; return) #write to log file about default setting
    inputDictionary["potential-unit"] == "hartree"  && (potential.potentialUnit = u"hartree"    ; return)
    inputDictionary["potential-unit"] == "ev"       && (potential.potentialUnit = u"eV"         ; return)
    inputDictionary["potential-unit"] == "kj/mol"   && (potential.potentialUnit = u"kJpermol"   ; return)
    inputDictionary["potential-unit"] == "kcal/mol" && (potential.potentialUnit = u"kcalpermol" ; return)

    @error "\nThe given potential-unit $(inputDictionary["potential-unit"]) was not recognised!\n" *
           "Valid opttions are:                                                                \n" * 
           "    - hartree                                                                      \n" *
           "    - ev                                                                           \n" *
           "    - kj/mol                                                                       \n" *
           "    - kcal/mol                                                                     \n"

    exit()
end

function checkCoordsUnits(potential::Potential)
    isempty(inputDictionary["coord-unit"])      && (potential.coordsUnit = u"angstrom"; return) #write to log file about default setting
    inputDictionary["coord-unit"] == "angstrom" && (potential.coordsUnit = u"angstrom"; return)
    inputDictionary["coord-unit"] == "nm"       && (potential.coordsUnit = u"nm"      ; return)
    inputDictionary["coord-unit"] == "bohr"     && (potential.coordsUnit = u"bohr"    ; return)

    @error "\nThe given coord-unit $(inputDictionary["coord-unit"]) was not recognised!\n" *
           "Valid opttions are:                                                        \n" * 
           "    - angstrom                                                             \n" *
           "    - nm                                                                   \n" *
           "    - bohr                                                                 \n"

    exit()
end

function checkMassUnits(potential::Potential)
    isempty(inputDictionary["mass-unit"])   && (potential.massUnit = u"u" ; return) #write to log file about default setting
    inputDictionary["mass-unit"] == "unit"  && (potential.massUnit = u"u" ; return)
    inputDictionary["mass-unit"] == "g/mol" && (potential.massUnit = u"u" ; return)
    inputDictionary["mass-unit"] == "me"    && (potential.massUnit = u"me"; return)

    @error "\nThe given mass-unit $(inputDictionary["mass-unit"]) was not recognised!\n" *
           "Valid opttions are:                                                      \n" * 
           "    - unit                                                               \n" *
           "    - g/mol                                                              \n" *
           "    - me            #stands for electron mass                            \n"

    exit()
end

function checkMass(potential::Potential)
    isempty(inputDictionary["reduced-mass"]) && (potential.mass = 1.0; return) #write to log file about default setting
    potential.mass = parse(Float64, inputDictionary["reduced-mass"])
end

function checkShiftPotential(potential::Potential)
    isempty(inputDictionary["shift-potential"]) && (potential.shift = false; return) #write to log file about default setting
end

function checkKPoints(potential::Potential)
    isempty(inputDictionary["k-points"]) && (potential.n_kpoints = -1; return) #write to log file about default setting
    potential.n_kpoints = parse(Int64, inputDictionary["k-points"])
    potential.n_kpoints < 2 && (@error "number of k-points has to be > 1"; exit())
end

function checkNDatapoints(potential::Potential)
    isempty(inputDictionary["datapoints"]) && (potential.n_datapoints = Vector(); return) #write to log file about default setting
    potential.n_datapoints = parse.(Int64, split(join(split(inputDictionary["datapoints"], ","), " ")))
    length(potential.n_datapoints) == 0 && (@error "n_datapoints not correctly defined!"; exit())
end

function checkBandStructure(potential::Potential)
    isempty(inputDictionary["band-structure"]) && (potential.bandStructure = false; return) #write to log file about default setting
    inputDictionary["band-structure"] ∈ ["on","true"] && (potential.bandStructure = true; return)
    @error ""; exit()
end

function checkStencil(system::System)
    isempty(inputDictionary["stencil"]) && (system.stencil = 9; return) #write to log file about default setting
    system.stencil = parse(Int64, inputDictionary["stencil"])

    system.stencil ∉ [3,5,7,9,11,13] && (@error "stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""; exit())
end

function checkPeriodicity(system::System)
    isempty(inputDictionary["periodic"]) && (system.periodic = false; return) #write to log file about default setting
    system.periodic = parse(Bool, inputDictionary["periodic"])
end

function checkNEigenvalues(output::Output)
    isempty(inputDictionary["n-eigenvalues"]) && (output.n_eigenvalues = 5; return) #write to log file about default setting
    output.n_eigenvalues = parse(Int64, inputDictionary["n-eigenvalues"])
    output.n_eigenvalues < 1 && (@error "number of eigenvalues has to be >= 1"; exit())
end