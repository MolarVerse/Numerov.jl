function checkInput(files::Files)
    checkPotentialFile(files)
    checkLogFileName(files)
    checkTimingsFileName(files)
end

function checkPotentialFile(files::Files)
    isempty(inputDictionary["potential-file"]) && (@error "No potential input file name given"; exit())
    files.PotentialFileName = inputDictionary["potential-file"]
end

function checkLogFileName(files::Files)
    files.logFileName = isempty(inputDictionary["output-file"]) ? "Numerov.out" : inputDictionary["output-file"]
    files.logFile     = open(files.logFileName, "w")
end

function checkTimingsFileName(files::Files)
    files.timingsFileName = isempty(inputDictionary["timings-file"]) ? "timings.out" : inputDictionary["timings-file"]
    files.timingsFile     = open(files.timingsFileName, "w")
end