module MyUnits

using Unitful
using PhysicalConstants

import PhysicalConstants.CODATA2018: h, ħ, N_A, c_0

@unit Jpermol   "J/mol"   myjoulepermol (1/ustrip(N_A))u"J" true
@unit calpermol "cal/mol" mycalpermol   (1/ustrip(N_A))u"cal" true
@unit me        "m_e"     electronMass   5.48579909065e-4*u"u" true

function __init__()
    # Guard against double registration: this is called explicitly at include
    # time (so that u"..." macros in later source files can resolve these units
    # during precompilation) and again automatically when the module is loaded.
    MyUnits in Unitful.unitmodules || Unitful.register(MyUnits)
end

end