function testcase(id)
    if id == "A"
        sc = basicSolarSail()                                       # Spacecraft struct
        X0 = [9222.7; 0.20; 0.573*pi/180; 0.00; 2.02; 0.0]          # Initial conditions
        XT = [26500.0, 0.75, 0.01*pi/180, 0.0*pi/180, 90*pi/180]    # Targets
        oetols = [10, 0.001, 0.01, 0.01, 0.01]                      # Tols
        Woe = [1.0, 1.0, 1.0, 0.0, 0.0]                             # Weights
        qlaw_type = Oguri                                           # Parameterization Type
    elseif id == "B"
        sc = basicSolarSail()
        X0 = [9222.7; 0.20; 0.573*pi/180; 0.00; 2.02; 0.0]
        XT = [26500.0, 0.10, 10.0*pi/180, 270.0*pi/180, 90.0*pi/180]
        oetols = [10, 0.001, 0.01, 0.01, 0.01]
        Woe = [22.0, 1.0, 1.0, 0.0, 0.0]
        qlaw_type = Oguri 
    elseif id == "C"
        sc = basicSolarSail()
        X0 = [26500; 0.70; 0.573*pi/180; 0.000; -4.179809375168203; 0.0]
        XT = [26700, 0.75, 0.2*pi/180, 30*pi/180, 90.0*pi/180] 
        oetols = [10, 0.001, 0.01, 0.01, 0.01]
        Woe = [1.0, 1.0, 1.0, 1.0, 0.0]
        qlaw_type = Oguri
    elseif id == "D"
        sc = basicSolarSail()
        X0 = [9222.7; 0.20; 0.573*pi/180; 0.00; 1.952; 0.0]
        XT = [26500.0, 0.75, 0.01*pi/180, 0.0*pi/180, 90*pi/180]
        oetols = [10, 0.001, 0.01, 0.01, 0.01]
        Woe = [1.0, 1.0, 1.0, 0.0, 0.0]
        qlaw_type = Keplerian
    elseif id == "E"
        sc = basicSolarSail()
        X0 = [26500; 0.70; 0.573*pi/180; 0.0; pi/4; 0.0]
        XT = [26700, 0.75, 0.2*pi/180, 30*pi/180, 90.0*pi/180] 
        oetols = [10, 0.001, 0.01, 0.01, 0.01]
        Woe = [1.0, 1.0, 1.0, 1.0, 0.0]
        qlaw_type = Keplerian
    elseif id == "F"
        sc = basicSolarSail()
        X0 = [26500; 0.210; 60.573*pi/180; 0.0; 0.0; 0.0]
        XT = [26700, 0.205, 60.0*pi/180, 30*pi/180, 30*pi/180]
        oetols = [10, 0.001, 0.01, 0.01, 0.01]
        Woe = [10.0, 10.0, 1.0, 1.0, 5.0]
        qlaw_type = Keplerian
    elseif id == "G"
        sc = basicSolarSail()
        X0 = [9222.7; 0.2; 0.573*pi/180; 0.0; 0.0; 0.0]
        XT = [26700, 0.7, 10.0*pi/180, 30*pi/180, 30*pi/180] 
        oetols = [10, 0.001, 0.01, 0.01, 0.01]
        Woe = [10.0, 10.0, 1.0, 1.0, 5.0]
        qlaw_type = Keplerian
    else
        error("Unrecognized case id, try again.")
    end
    
    return sc, X0, XT, oetols, Woe, qlaw_type
end