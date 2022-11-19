function buildLaplace(system::System1D)

    n_datapoints = system.n_datapoints
    n_datapoints < system.stencil && (@error "The number of datapoints has at least to be equal to the stencil size!"; exit())

    if system.stencil == 3

        system.laplace = spdiagm(
            -1 => ones(n_datapoints-1),
             0 => ones(n_datapoints  )*(-2),
             1 => ones(n_datapoints-1)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-1 => ones(1),
                -n_datapoints+1 => ones(1)
            )

        end

    elseif system.stencil == 5
        
        system.laplace = spdiagm(
            -2 => ones(n_datapoints-2)*(-1/12),
            -1 => ones(n_datapoints-1)*( 4/3 ),
             0 => ones(n_datapoints  )*(-5/2 ),
             1 => ones(n_datapoints-1)*( 4/3 ),
             2 => ones(n_datapoints-2)*(-1/12)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-2 => ones(2)*(-1/12),
                 n_datapoints-1 => ones(1)*( 4/3),
                -n_datapoints+1 => ones(1)*( 4/3),
                -n_datapoints+2 => ones(2)*(-1/12)
            )
            
        end

    elseif system.stencil == 7
        
        system.laplace = spdiagm(
            -3 => ones(n_datapoints-3)*(  1/90),
            -2 => ones(n_datapoints-2)*( -3/20),
            -1 => ones(n_datapoints-1)*(  3/2 ),
             0 => ones(n_datapoints  )*(-49/18),
             1 => ones(n_datapoints-1)*(  3/2 ),
             2 => ones(n_datapoints-2)*( -3/20),
             3 => ones(n_datapoints-3)*(  1/90)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-3 => ones(3)*(  1/90),
                 n_datapoints-2 => ones(2)*( -3/20),
                 n_datapoints-1 => ones(1)*(  3/2 ),
                -n_datapoints+1 => ones(1)*(  3/2 ),
                -n_datapoints+2 => ones(2)*( -3/20),
                -n_datapoints+3 => ones(3)*(  1/90)
            )
            
        end

    elseif system.stencil == 9
        
        system.laplace = spdiagm(
            -4 => ones(n_datapoints-4)*(  -1/560),
            -3 => ones(n_datapoints-3)*(   8/315),
            -2 => ones(n_datapoints-2)*(  -1/5  ),
            -1 => ones(n_datapoints-1)*(   8/5  ),
             0 => ones(n_datapoints  )*(-205/72 ),
             1 => ones(n_datapoints-1)*(   8/5  ),
             2 => ones(n_datapoints-2)*(  -1/5  ),
             3 => ones(n_datapoints-3)*(   8/315),
             4 => ones(n_datapoints-4)*(  -1/560)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-4 => ones(4)*(  -1/560),
                 n_datapoints-3 => ones(3)*(   8/315),
                 n_datapoints-2 => ones(2)*(  -1/5  ),
                 n_datapoints-1 => ones(1)*(   8/5  ),
                -n_datapoints+1 => ones(1)*(   8/5  ),
                -n_datapoints+2 => ones(2)*(  -1/5  ),
                -n_datapoints+3 => ones(3)*(   8/315),
                -n_datapoints+4 => ones(4)*(  -1/560)
            )
            
        end

    elseif system.stencil == 11
        
        system.laplace = spdiagm(
            -5 => ones(n_datapoints-5)*(     8/25200),
            -4 => ones(n_datapoints-4)*(  -125/25200),
            -3 => ones(n_datapoints-3)*(  1000/25200),
            -2 => ones(n_datapoints-2)*( -6000/25200),
            -1 => ones(n_datapoints-1)*( 42000/25200),
             0 => ones(n_datapoints  )*(-73766/25200),
             1 => ones(n_datapoints-1)*( 42000/25200),
             2 => ones(n_datapoints-2)*( -6000/25200),
             3 => ones(n_datapoints-3)*(  1000/25200),
             4 => ones(n_datapoints-4)*(  -125/25200),
             5 => ones(n_datapoints-5)*(     8/25200)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-5 => ones(5)*(     8/25200),
                 n_datapoints-4 => ones(4)*(  -125/25200),
                 n_datapoints-3 => ones(3)*(  1000/25200),
                 n_datapoints-2 => ones(2)*( -6000/25200),
                 n_datapoints-1 => ones(1)*( 42000/25200),
                -n_datapoints+1 => ones(1)*( 42000/25200),
                -n_datapoints+2 => ones(2)*( -6000/25200),
                -n_datapoints+3 => ones(3)*(  1000/25200),
                -n_datapoints+4 => ones(4)*(  -125/25200),
                -n_datapoints+5 => ones(5)*(     8/25200)
            )
            
        end

    elseif system.stencil == 13

        system.laplace = spdiagm(
            -6 => ones(n_datapoints-6)*(     -50/831600),
            -5 => ones(n_datapoints-5)*(     864/831600),
            -4 => ones(n_datapoints-4)*(   -7425/831600),
            -3 => ones(n_datapoints-3)*(   44000/831600),
            -2 => ones(n_datapoints-2)*( -222750/831600),
            -1 => ones(n_datapoints-1)*( 1425600/831600),
             0 => ones(n_datapoints  )*(-2480478/831600),
             1 => ones(n_datapoints-1)*( 1425600/831600),
             2 => ones(n_datapoints-2)*( -222750/831600),
             3 => ones(n_datapoints-3)*(   44000/831600),
             4 => ones(n_datapoints-4)*(   -7425/831600),
             5 => ones(n_datapoints-5)*(     864/831600),
             6 => ones(n_datapoints-6)*(     -50/831600)
        )

        if system.periodic

            system.laplace += spdiagm(
                 n_datapoints-6 => ones(6)*(     -50/831600),
                 n_datapoints-5 => ones(5)*(     864/831600),
                 n_datapoints-4 => ones(4)*(   -7425/831600),
                 n_datapoints-3 => ones(3)*(   44000/831600),
                 n_datapoints-2 => ones(2)*( -222750/831600),
                 n_datapoints-1 => ones(1)*( 1425600/831600),
                -n_datapoints+1 => ones(1)*( 1425600/831600),
                -n_datapoints+2 => ones(2)*( -222750/831600),
                -n_datapoints+3 => ones(3)*(   44000/831600),
                -n_datapoints+4 => ones(4)*(   -7425/831600),
                -n_datapoints+5 => ones(5)*(     864/831600),
                -n_datapoints+6 => ones(6)*(     -50/831600)
            )
            
        end

    end

end

function buildLaplace(system::System2D)

    n_datapoints   = system.n_datapoints
    system.laplace = zeros(n_datapoints[1]*n_datapoints[2], n_datapoints[1]*n_datapoints[2])

    n_datapoints[1] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())
    n_datapoints[2] < system.stencil && (@error "The number of datapoints in each dimension has at least to be equal to the stencil size!"; exit())

    stencil = zeros(system.stencil, system.stencil)

    if system.stencil == 3

        stencil[:,1] = [1.0,  0.0, 1.0]
        stencil[:,2] = [0.0, -4.0, 0.0]
        stencil[:,3] = [1.0,  0.0, 1.0]

    elseif system.stencil == 5
        
        stencil[:,1] = [ 0.0,  0.0,  -1.0,  0.0,  0.0]
        stencil[:,2] = [ 0.0,  0.0,  16.0,  0.0,  0.0]
        stencil[:,3] = [-1.0, 16.0, -60.0, 16.0, -1.0]
        stencil[:,4] = [ 0.0,  0.0,  16.0,  0.0,  0.0]
        stencil[:,5] = [ 0.0,  0.0,  -1.0,  0.0,  0.0]

        stencil /= 6.0

    elseif system.stencil == 7
        
        stencil[:,1] = [ 0.0,    0.0,    0.0,    16.0,    0.0,    0.0,  0.0]
        stencil[:,2] = [ 0.0,   -5.0,   20.0,  -246.0,   20.0,   -5.0,  0.0]
        stencil[:,3] = [ 0.0,   20.0,  -80.0,  2280.0,  -80.0,   20.0,  0.0]
        stencil[:,4] = [16.0, -246.0, 2280.0, -8020.0, 2280.0, -246.0, 16.0]
        stencil[:,5] = [ 0.0,   20.0,  -80.0,  2280.0,  -80.0,   20.0,  0.0]
        stencil[:,6] = [ 0.0,   -5.0,   20.0,  -246.0,   20.0,   -5.0,  0.0]
        stencil[:,7] = [ 0.0,    0.0,    0.0,    16.0,    0.0,    0.0,  0.0]

        stencil /= 720.0

    elseif system.stencil == 9
        
        stencil[:,1] = [    0.0,     0.0,       0.0,       0.0,    -1620.0,       0.0,       0.0,     0.0,     0.0]
        stencil[:,2] = [    0.0,   -56.0,     231.0,    -420.0,    23530.0,    -420.0,     231.0,   -56.0,     0.0]
        stencil[:,3] = [    0.0,   231.0,    -756.0,     945.0,  -182280.0,     945.0,    -756.0,   231.0,     0.0]
        stencil[:,4] = [    0.0,  -420.0,     945.0,       0.0,  1450470.0,       0.0,     945.0,  -420.0,     0.0]
        stencil[:,5] = [-1620.0, 23530.0, -182280.0, 1450470.0, -5163200.0, 1450470.0, -182280.0, 23530.0, -1620.0]
        stencil[:,6] = [    0.0,  -420.0,     945.0,       0.0,  1450470.0,       0.0,     945.0,  -420.0,     0.0]
        stencil[:,7] = [    0.0,   231.0,    -756.0,     945.0,  -182280.0,     945.0,    -756.0,   231.0,     0.0]
        stencil[:,8] = [    0.0,   -56.0,     231.0,    -420.0,    23530.0,    -420.0,     231.0,   -56.0,     0.0]
        stencil[:,9] = [    0.0,     0.0,       0.0,       0.0,    -1620.0,       0.0,       0.0,     0.0,     0.0]

        stencil /= 453600.0

    elseif system.stencil == 11
        
        stencil[:, 1] = [    0.0,       0.0,       0.0,         0.0,        0.0,      16128.0,        0.0,         0.0,       0.0,       0.0,     0.0]
        stencil[:, 2] = [    0.0,     -81.0,     466.0,     -1281.0,     2226.0,    -254660.0,     2226.0,     -1281.0,     466.0,     -81.0,     0.0]
        stencil[:, 3] = [    0.0,     466.0,   -2468.0,      6328.0,   -10556.0,    2028460.0,   -10556.0,      6328.0,   -2468.0,     466.0,     0.0]
        stencil[:, 4] = [    0.0,   -1281.0,    6328.0,    -15288.0,    24696.0,  -12124910.0,    24696.0,    -15288.0,    6328.0,   -1281.0,     0.0]
        stencil[:, 5] = [    0.0,    2226.0,  -10556.0,     24696.0,   -39396.0,   84718060.0,   -39396.0,     24696.0,  -10556.0,    2226.0,     0.0]
        stencil[:, 6] = [16128.0, -254660.0, 2028460.0, -12124910.0, 84718060.0, -297478412.0, 84718060.0, -12124910.0, 2028460.0, -254660.0, 16128.0]
        stencil[:, 7] = [    0.0,    2226.0,  -10556.0,     24696.0,   -39396.0,   84718060.0,   -39396.0,     24696.0,  -10556.0,    2226.0,     0.0]
        stencil[:, 8] = [    0.0,   -1281.0,    6328.0,    -15288.0,    24696.0,  -12124910.0,    24696.0,    -15288.0,    6328.0,   -1281.0,     0.0]
        stencil[:, 9] = [    0.0,     466.0,   -2468.0,      6328.0,   -10556.0,    2028460.0,   -10556.0,      6328.0,   -2468.0,     466.0,     0.0]
        stencil[:,10] = [    0.0,     -81.0,     466.0,     -1281.0,     2226.0,    -254660.0,     2226.0,     -1281.0,     466.0,     -81.0,     0.0]
        stencil[:,11] = [    0.0,       0.0,       0.0,         0.0,        0.0,      16128.0,        0.0,         0.0,       0.0,       0.0,     0.0]
        stencil[:,10] = [    0.0,     -81.0,     466.0,     -1281.0,     2226.0,    -254660.0,     2226.0,     -1281.0,     466.0,     -81.0,     0.0]
        stencil[:, 9] = [    0.0,     466.0,   -2468.0,      6328.0,   -10556.0,    2028460.0,   -10556.0,      6328.0,   -2468.0,     466.0,     0.0]
        stencil[:, 8] = [    0.0,   -1281.0,    6328.0,    -15288.0,    24696.0,  -12124910.0,    24696.0,    -15288.0,    6328.0,   -1281.0,     0.0]
        stencil[:, 7] = [    0.0,    2226.0,  -10556.0,     24696.0,   -39396.0,   84718060.0,   -39396.0,     24696.0,  -10556.0,    2226.0,     0.0]
        stencil[:, 6] = [16128.0, -254660.0, 2028460.0, -12124910.0, 84718060.0, -297478412.0, 84718060.0, -12124910.0, 2028460.0, -254660.0, 16128.0]
        stencil[:, 5] = [    0.0,    2226.0,  -10556.0,     24696.0,   -39396.0,   84718060.0,   -39396.0,     24696.0,  -10556.0,    2226.0,     0.0]
        stencil[:, 4] = [    0.0,   -1281.0,    6328.0,    -15288.0,    24696.0,  -12124910.0,    24696.0,    -15288.0,    6328.0,   -1281.0,     0.0]
        stencil[:, 3] = [    0.0,     466.0,   -2468.0,      6328.0,   -10556.0,    2028460.0,   -10556.0,      6328.0,   -2468.0,     466.0,     0.0]
        stencil[:, 2] = [    0.0,     -81.0,     466.0,     -1281.0,     2226.0,    -254660.0,     2226.0,     -1281.0,     466.0,     -81.0,     0.0]
        stencil[:, 1] = [    0.0,       0.0,       0.0,         0.0,        0.0,      16128.0,        0.0,         0.0,       0.0,       0.0,     0.0]

        stencil /= 25401600.0

    elseif system.stencil == 13

        stencil[:, 1] = [0.0, 0.0,0.0,0.0,0.0,0.0,-58517879280.0,0.0,0.0,0.0,0.0,0.0,0.0]
        stencil[:, 2] = [0.0,-49045709.0,373590360.0,-1358549104.0,3109421304.0,-4962046320.0,1016963667642.0,-4962046320.0,3109421304.0,-1358549104.0,373590360.0,-49045709.0,0.0]
        stencil[:, 3] = [0.0,373590360.0,-2769597765.0,9878325160.0,-22333998060.0,35404871040.0,-8731046409470.0,35404871040.0,-22333998060.0,9878325160.0,-2769597765.0,373590360.0,0.0]
        stencil[:, 4] = [0.0,-1358549104.0,9878325160.0,-34830116640.0,78335831640.0,-123961751760.0,51640086176448.0,-123961751760.0,78335831640.0,-34830116640.0,9878325160.0,-1358549104.0,0.0]
        stencil[:, 5] = zeros(13)

        # stencil[ 52] =                 0.0,
        # stencil[ 53] =        3109421304.0/486656045568000.0,
        # stencil[ 54] =      -22333998060.0/486656045568000.0,
        # stencil[ 55] =       78335831640.0/486656045568000.0,
        # stencil[ 56] =     -176123906640.0/486656045568000.0,
        # stencil[ 57] =      279081578160.0/486656045568000.0,
        # stencil[ 58] =  -261377128691448.0/486656045568000.0,
        # stencil[ 59] =      279081578160.0/486656045568000.0,
        # stencil[ 60] =     -176123906640.0/486656045568000.0,
        # stencil[ 61] =       78335831640.0/486656045568000.0,
        # stencil[ 62] =      -22333998060.0/486656045568000.0,
        # stencil[ 63] =        3109421304.0/486656045568000.0,
        # stencil[ 64] =                 0.0,

        # stencil[ 65] =                 0.0,
        # stencil[ 66] =       -4962046320.0/486656045568000.0,
        # stencil[ 67] =       35404871040.0/486656045568000.0,
        # stencil[ 68] =     -123961751760.0/486656045568000.0,
        # stencil[ 69] =      279081578160.0/486656045568000.0,
        # stencil[ 70] =     -443097325440.0/486656045568000.0,
        # stencil[ 71] =  1669039023109920.0/486656045568000.0,
        # stencil[ 72] =     -443097325440.0/486656045568000.0,
        # stencil[ 73] =      279081578160.0/486656045568000.0,
        # stencil[ 74] =     -123961751760.0/486656045568000.0,
        # stencil[ 75] =       35404871040.0/486656045568000.0,
        # stencil[ 76] =       -4962046320.0/486656045568000.0,
        # stencil[ 77] =                 0.0,

        # stencil[ 78] =      -58517879280.0/486656045568000.0,
        # stencil[ 79] =     1016963667642.0/486656045568000.0,
        # stencil[ 80] =    -8731046409470.0/486656045568000.0,
        # stencil[ 81] =    51640086176448.0/486656045568000.0,
        # stencil[ 82] =  -261377128691448.0/486656045568000.0,
        # stencil[ 83] =  1669039023109920.0/486656045568000.0,
        # stencil[ 84] = -5806918679824232.0/486656045568000.0,
        # stencil[ 85] =  1669039023109920.0/486656045568000.0,
        # stencil[ 86] =  -261027103261848.0/486656045568000.0,
        # stencil[ 87] =    51640086176448.0/486656045568000.0,
        # stencil[ 88] =    -8731046409470.0/486656045568000.0,
        # stencil[ 89] =     1016963667642.0/486656045568000.0,
        # stencil[ 90] =      -58517879280.0/486656045568000.0,

        # stencil[ 91] =                 0.0,
        # stencil[ 92] =       -4962046320.0/486656045568000.0,
        # stencil[ 93] =       35404871040.0/486656045568000.0,
        # stencil[ 94] =     -123961751760.0/486656045568000.0,
        # stencil[ 95] =      279081578160.0/486656045568000.0,
        # stencil[ 96] =     -443097325440.0/486656045568000.0,
        # stencil[ 97] =  1669039023109920.0/486656045568000.0,
        # stencil[ 98] =     -443097325440.0/486656045568000.0,
        # stencil[ 99] =      279081578160.0/486656045568000.0,
        # stencil[100] =     -123961751760.0/486656045568000.0,
        # stencil[101] =       35404871040.0/486656045568000.0,
        # stencil[102] =       -4962046320.0/486656045568000.0,
        # stencil[103] =                 0.0,

        # stencil[104] =                 0.0,
        # stencil[105] =        3109421304.0/486656045568000.0,
        # stencil[106] =      -22333998060.0/486656045568000.0,
        # stencil[107] =       78335831640.0/486656045568000.0,
        # stencil[108] =     -176123906640.0/486656045568000.0,
        # stencil[109] =      279081578160.0/486656045568000.0,
        # stencil[110] =  -261377128691448.0/486656045568000.0,
        # stencil[111] =      279081578160.0/486656045568000.0,
        # stencil[112] =     -176123906640.0/486656045568000.0,
        # stencil[113] =       78335831640.0/486656045568000.0,
        # stencil[114] =      -22333998060.0/486656045568000.0,
        # stencil[115] =        3109421304.0/486656045568000.0,
        # stencil[116] =                 0.0,

        # stencil[117] =                 0.0,
        # stencil[118] =       -1358549104.0/486656045568000.0,
        # stencil[119] =        9878325160.0/486656045568000.0,
        # stencil[120] =      -34830116640.0/486656045568000.0,
        # stencil[121] =       78335831640.0/486656045568000.0,
        # stencil[122] =     -123961751760.0/486656045568000.0,
        # stencil[123] =    51640086176448.0/486656045568000.0,
        # stencil[124] =     -123961751760.0/486656045568000.0,
        # stencil[125] =       78335831640.0/486656045568000.0,
        # stencil[126] =      -34830116640.0/486656045568000.0,
        # stencil[127] =        9878325160.0/486656045568000.0,
        # stencil[128] =       -1358549104.0/486656045568000.0,
        # stencil[129] =                 0.0,

        # stencil[130] =                 0.0,
        # stencil[131] =         373590360.0/486656045568000.0,
        # stencil[132] =       -2769597765.0/486656045568000.0,
        # stencil[133] =        9878325160.0/486656045568000.0,
        # stencil[134] =      -22333998060.0/486656045568000.0,
        # stencil[135] =       35404871040.0/486656045568000.0,
        # stencil[136] =    -8731046409470.0/486656045568000.0,
        # stencil[137] =       35404871040.0/486656045568000.0,
        # stencil[138] =      -22333998060.0/486656045568000.0,
        # stencil[139] =        9878325160.0/486656045568000.0,
        # stencil[140] =       -2769597765.0/486656045568000.0,
        # stencil[141] =         373590360.0/486656045568000.0,
        # stencil[142] =                 0.0,

        # stencil[143] =                 0.0,
        # stencil[144] =         -49045709.0/486656045568000.0,
        # stencil[145] =         373590360.0/486656045568000.0,
        # stencil[146] =       -1358549104.0/486656045568000.0,
        # stencil[147] =        3109421304.0/486656045568000.0,
        # stencil[148] =       -4962046320.0/486656045568000.0,
        # stencil[149] =     1016963667642.0/486656045568000.0,
        # stencil[150] =       -4962046320.0/486656045568000.0,
        # stencil[151] =        3109421304.0/486656045568000.0,
        # stencil[152] =       -1358549104.0/486656045568000.0,
        # stencil[153] =         373590360.0/486656045568000.0,
        # stencil[154] =         -49045709.0/486656045568000.0,
        # stencil[155] =                 0.0,

        # stencil[156] =                 0.0,
        # stencil[157] =                 0.0,
        # stencil[158] =                 0.0,
        # stencil[159] =                 0.0,
        # stencil[160] =                 0.0,
        # stencil[161] =                 0.0,
        # stencil[162] =      -58517879280.0/486656045568000.0,
        # stencil[163] =                 0.0,
        # stencil[164] =                 0.0,
        # stencil[165] =                 0.0,
        # stencil[166] =                 0.0,
        # stencil[167] =                 0.0,
        # stencil[168] =                 0.0,

        stencil /= 486656045568000.0

    end

    for x_1 in 1:n_datapoints[1]
        for stencil_index in 1:system.stencil
    
            super_matrix_index = x_1 - system.stencil ÷ 2 - 1 + stencil_index
    
            if super_matrix_index < 1 && system.periodic == false
                continue
            elseif super_matrix_index < 1
                super_matrix_index += n_datapoints[1]
            end
    
            matrix = zeros(n_datapoints[2], n_datapoints[2])
    
            for j in 1:system.stencil
    
                sub_matrix_index = j - system.stencil ÷ 2 - 1
    
                if sub_matrix_index == 0
                    matrix += spdiagm(sub_matrix_index => ones(n_datapoints[2] - sub_matrix_index) * stencil[stencil_index, j])
                else
                    matrix += spdiagm(sub_matrix_index => ones(n_datapoints[2] - abs(sub_matrix_index)) * stencil[stencil_index, j])
                end
    
                if system.periodic && sub_matrix_index != 0
                    matrix += spdiagm(sign(sub_matrix_index)*n_datapoints[2] - sub_matrix_index => ones(abs(sub_matrix_index)) * stencil[stencil_index, j])
                end
            end
    
            i_1 = 1 + (x_1 - 1) * n_datapoints[2]
            j_1 = x_1 * n_datapoints[2]
    
            i_2 = 1 + (super_matrix_index - 1) * n_datapoints[2]
            j_2 = (super_matrix_index) * n_datapoints[2]
    
            if(i_2 > n_datapoints[2]*n_datapoints[1])
                continue
            end
    
            system.laplace[i_1:j_1, i_2:j_2] += matrix #add system.laplace
    
        end
    end
    
    for i in 1:n_datapoints[1]*n_datapoints[2]
        for j in i:n_datapoints[1]*n_datapoints[2]
            system.laplace[j,i] = system.laplace[i,j]
        end
    end
end