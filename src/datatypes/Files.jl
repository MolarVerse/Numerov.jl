mutable struct Files

    eigenvalueFileName             ::String
    eigenvectorFileName            ::String
    eigenvectorShiftedFileName     ::String
    imag_eigenvectorFileName       ::String
    imag_eigenvectorShiftedFileName::String
    frequencyFileName              ::String
    bandStructureFileName          ::String

    inputFileName    ::String
    PotentialFileName::String

    logFileName::String
    logFile    ::IOStream

    timingsFileName::String
    timingsFile    ::IOStream

    k_pointsFileName::String

    to::TimerOutput

    Files() = new()
end