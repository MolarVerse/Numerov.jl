function checkInput(system::System)
    checkStencil(system)
    checkStencilLaplace(system)
    checkStencilNabla(system)
    checkSolver(system)
end

function checkStencil(system::System)
    isempty(inputDictionary["stencil"]) && (system.stencil = 9; return) #write to log file about default setting
    system.stencil = parse(Int64, inputDictionary["stencil"])
    system.stencil ∉ [3,5,7,9,11,13] && throw(ArgumentError("stencil $(inputDictionary["stencil"]) not available -- possible entries are \"3,5,7,9,11,13\""))
end

function checkStencilLaplace(system::System)
    isempty(inputDictionary["stencil-laplace"]) && (system.stencilΔ = 0; return) #write to log file about default setting
    system.stencilΔ = parse(Int64, inputDictionary["stencil-laplace"])
    system.stencilΔ ∉ [3,5,7,9,11,13] && throw(ArgumentError("stencil-laplace $(inputDictionary["stencil-laplace"]) not available -- possible entries are \"3,5,7,9,11,13\""))
end

function checkStencilNabla(system::System)
    isempty(inputDictionary["stencil-nabla"]) && (system.stencil∇ = 0; return) #write to log file about default setting
    system.stencil∇ = parse(Int64, inputDictionary["stencil-nabla"])
    system.stencil∇ ∉ [3,5,7,9,11] && throw(ArgumentError("stencil-nabla $(inputDictionary["stencil-nabla"]) not available -- possible entries are \"3,5,7,9,11\""))
end

function checkSolver(system::System)
    solver = lowercase(inputDictionary["solver"])
    isempty(solver)    && (system.solver = ARPACK ; return) #write to log file about default setting
    solver == "arpack" && (system.solver = ARPACK ; return)
    solver == "krylov" && (system.solver = KRYLOV ; return)
    solver == "lobpcg" && (system.solver = LOBPCG; return)
    solver == "cuda"   && throw(ArgumentError("the cuda solver is not implemented -- use arpack, krylov or lu"))
    solver == "lu"     && (system.solver = LU     ; return)

    throw(ArgumentError("\nThe given solver $(inputDictionary["solver"]) was not recognised!\n" *
                        "Valid options are:                                               \n" *
                        "    - arpack                                                     \n" *
                        "    - krylov                                                     \n" *
                        "    - lobpcg                                                     \n" *
                        "    - lu                                                         \n"))
end
