mutable struct Files

    eigenvectorFileName            ::String
    eigenvectorShiftedFileName     ::String
    imag_eigenvectorFileName       ::String
    imag_eigenvectorShiftedFileName::String
    frequencyFileName              ::String

    logFileName    ::String
    timingsFileName::String

    logFile    ::IOStream
    timingsFile::IOStream

    to::TimerOutput

    Files() = new()
end