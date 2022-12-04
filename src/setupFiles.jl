function setupFiles(files::Files)
    checkLogFile(files)
end

function checkLogFile(files::Files)
    files.logFileName = isempty(inputDictionary["output-file"]) ? "Numerov.out" : inputDictionary["output-file"]
    files.logFile = open(files.logFileName, "w")
end