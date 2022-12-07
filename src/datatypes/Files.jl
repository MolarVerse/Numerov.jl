mutable struct Files
    logFileName    ::String
    timingsFileName::String

    logFile    ::IOStream
    timingsFile::IOStream

    to::TimerOutput

    Files() = new()
end