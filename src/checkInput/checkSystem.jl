function checkInput(system::System)
    checkStencil(system)
    checkStencilLaplace(system)
    checkStencilNabla(system)
    checkSolver(system)
end

function checkStencil(system::System)
    isempty(inputDictionary["stencil"]) && (system.stencil = 9; return) #write to log file about default setting
    system.stencil = parse(Int64, inputDictionary["stencil"])
    system.stencil ∉ [3,5,7,9,11,13] && (@error "stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""; exit())
end

function checkStencilLaplace(system::System)
    isempty(inputDictionary["stencil-laplace"]) && (system.stencilΔ = 0; return) #write to log file about default setting
    system.stencilΔ = parse(Int64, inputDictionary["stencil-laplace"])
    system.stencilΔ ∉ [3,5,7,9,11,13] && (@error "stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""; exit())
end

function checkStencilNabla(system::System)
    isempty(inputDictionary["stencil-nabla"]) && (system.stencil∇ = 0; return) #write to log file about default setting
    system.stencil∇ = parse(Int64, inputDictionary["stencil-nabla"])
    system.stencil∇ ∉ [3,5,7,9,11] && (@error "stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""; exit())
end

function checkSolver(system::System)
    isempty(inputDictionary["solver"])    && (system.solver = ARPACK ; return) #write to log file about default setting
    inputDictionary["solver"] == "arpack" && (system.solver = ARPACK ; return)
    inputDictionary["solver"] == "krylov" && (system.solver = KRYLOV ; return)
    inputDictionary["solver"] == "cuda"   && (system.solver = GPU    ; return)
    inputDictionary["solver"] == "lu"     && (system.solver = LU     ; return)

    @error "\nThe given solver $(inputDictionary["solver"]) was not recognised!\n" *
           "Valid opttions are:                                                      \n" * 
           "    - arpack                                                             \n" *
           "    - krylov                                                             \n" *
           "    - cuda                                                               \n" *
           "    - LU                                                                 \n"

    exit()
end