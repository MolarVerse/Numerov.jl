function checkInput(output::Output)
    checkNEigenvalues(output)
end

function checkNEigenvalues(output::Output)
    isempty(inputDictionary["n-eigenvalues"]) && (output.n_eigenvalues = 5; return) #write to log file about default setting
    output.n_eigenvalues = parse(Int64, inputDictionary["n-eigenvalues"])
    output.n_eigenvalues < 1 && throw(ArgumentError("number of eigenvalues has to be >= 1"))
end