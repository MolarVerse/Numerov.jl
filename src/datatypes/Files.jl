mutable struct Files
    logFileName::String

    logFile::IOStream

    Files() = new()
end