function checkInput(system::System)
    checkStencil(system)
    checkPeriodicity(system)
    checkSolver(system)
end

function checkStencil(system::System)
    isempty(inputDictionary["stencil"]) && (system.stencil = 9; return) #write to log file about default setting
    system.stencil = parse(Int64, inputDictionary["stencil"])

    system.stencil ∉ [3,5,7,9,11,13] && (@error "stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""; exit())
end

function checkPeriodicity(system::System)
    isempty(inputDictionary["periodic"]) && (system.periodic = [false]; return) #write to log file about default setting
    system.periodic = parse.(Bool, split(join(split(inputDictionary["periodic"], ","), " ")))
end

function checkSolver(system::System)
    isempty(inputDictionary["solver"])    && (system.solver = ARPACK ; return) #write to log file about default setting
    inputDictionary["solver"] == "arpack" && (system.solver = ARPACK ; return)
    inputDictionary["solver"] == "krylov" && (system.solver = KRYLOV ; return)
    inputDictionary["solver"] == "cuda"   && (system.solver = GPU    ; return)

    @error "\nThe given solver $(inputDictionary["solver"]) was not recognised!\n" *
           "Valid opttions are:                                                      \n" * 
           "    - arpack                                                             \n" *
           "    - krylov                                                             \n" *
           "    - cuda                                                               \n"

    exit()
end