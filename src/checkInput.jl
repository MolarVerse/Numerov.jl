function checkInput(potential::Potential)
    checkPotentialFile(potential)
    checkPotentialUnits(potential)
    checkCoordsUnits(potential)
end

function checkInput(system::System)
    checkStencil(system)
end

function checkPotentialFile(potential::Potential)
    isempty(inputDictionary["potential-file"]) && (@error "No potential input file name given"; exit())
    potential.file = inputDictionary["potential-file"]
end

function checkPotentialUnits(potential::Potential)
    isempty(inputDictionary["potential-unit"]) && (potential.potentialUnit = u"hartree"; return) #write to log file about default setting
end

function checkCoordsUnits(potential::Potential)
    isempty(inputDictionary["coord-unit"]) && (potential.coordsUnit = u"bohr"; return) #write to log file about default setting
end

function checkStencil(system::System)
    isempty(inputDictionary["stencil"]) && (system.stencil = 9; return) #write to log file about default setting
    system.stencil = parse(Int64, inputDictionary["stencil"])

    system.stencil ∉ [3,5,7,9,11,13] && (@error "stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""; exit())
end