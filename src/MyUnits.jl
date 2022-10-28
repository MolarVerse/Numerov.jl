module MyUnits

using Unitful
using PhysicalConstants

import PhysicalConstants.CODATA2018: ħ, N_A

@unit Jpermol   "J/mol"   myjoulepermol (1/ustrip(N_A))u"J" true;
@unit calpermol "cal/mol" mycalpermol   (1/ustrip(N_A))u"cal" true;

function __init__()
    Unitful.register(MyUnits)
end

end