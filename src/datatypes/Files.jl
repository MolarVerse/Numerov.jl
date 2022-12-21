mutable struct Files

    eigenvectorFileName            ::String
    eigenvectorShiftedFileName     ::String
    imag_eigenvectorFileName       ::String
    imag_eigenvectorShiftedFileName::String
    frequencyFileName              ::String

    inputFileName    ::String
    PotentialFileName::String

    logFileName::String
    logFile    ::IOStream

    timingsFileName::String
    timingsFile    ::IOStream

    to::TimerOutput

    Files() = new()
end