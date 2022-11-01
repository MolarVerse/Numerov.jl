__precompile__(true)
module MyUnits

using Unitful
using PhysicalConstants

import PhysicalConstants.CODATA2018: h, ħ, N_A, c_0

@unit Jpermol   "J/mol"   myjoulepermol (1/ustrip(N_A))u"J" true
@unit calpermol "cal/mol" mycalpermol   (1/ustrip(N_A))u"cal" true
@unit me        "m_e"     electronMass   5.48579909065e-4*u"u" true

const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(MyUnits)
end

end