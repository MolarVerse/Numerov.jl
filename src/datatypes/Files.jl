mutable struct Files
    logFileName    ::String
    timingsFileName::String

    logFile    ::IOStream
    timingsFile::IOStream

    Files() = new()
end