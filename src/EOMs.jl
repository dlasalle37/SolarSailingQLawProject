function solarSailEOM_cartesian!(dx, x, p, t)
    #Unpack Parameters and current step values
        # paramter set p is organized as [mu, ::basicSolarSail, ::TwoBodyEphemeride]
    mu = p[1]
    sc = p[2]
    eph = p[3]
    epoch = eph.t0
    rVec = x[1:3]
    vVec = x[4:6]
    r = norm(rVec)

    # Compute SRP acceleration
    # A constant SRP accel for now (for fun!) (export all of this to a calculation function LATER)
    C1 = sc.C[1]; C2 = sc.C[2]; C3=sc.C[3];
    A = sc.areaParam
    nHat = [1; 0; 0]  # sail vector aligned with inertial x
    G0 = 1.02E14 # [kgkm/s^2]
    sunDist = 1.46E8
    (sHat, sunDist) = get_sunlight_direction(epoch+t)  # sHat expressed in ECI frame
    α = pi - acos(dot(nHat, sHat))  # alpha is the angle b/w the anti-sunlight direction and nHat
    a_SRP = A*G0/(sunDist^2)*cos(α)*(-(C1*cos(α)+C2)*nHat + C3*sHat)

    # Only thrust in velocity direction
    if dot(a_SRP, vVec) <= 0
        α = pi/2
        a_SRP = A*G0/(sunDist^2)*cos(α)*(-(C1*cos(α)+C2)*nHat + C3*sHat)
    end

    # EOM
    dx[1:3] .= vVec
    dx[4] = -mu/r^3 * rVec[1] + a_SRP[1]
    dx[5] = -mu/r^3 * rVec[2] + a_SRP[2]
    dx[6] = -mu/r^3 * rVec[3] + a_SRP[3]
end

"""
aSRP: calculating SRP acceleration in the Hill (Sun centered) frame
Notes:
    - SRP acceleration vector is calculated in the sun-centered hill frame
inputs:
    u: control variables [alpha, beta]
    sc: basicSolarSail
    eph: ephemeride object for simulation
    trueAnom_Earth: earth current true anomaly

outputs:
    aSRP: SRP acceleration vector in Hill (sun-centered) frame
"""
function aSRP(u, sc::basicSolarSail, eph::TwoBodyEphemeride, trueAnom_Earth)
    # Unpack Inputs:
    α = u[1]
    β = u[2]
    C1 = sc.C[1]
    C2 = sc.C[2]
    C3 = sc.C[3]
    d = distance_to_sun(eph, trueAnom_Earth)
    C1 = sc.C[1]; C2 = sc.C[2]; C3 = sc.C[3]
    G0 = 1.02E14 # solar flux constant at earth [kgkm/s^2]
    a1 = C1*(cos(α))^2+C2*cos(α)+C3  # sun frame r-direction, not scaled
    a2 = -(C1*cos(α)+C2)*sin(α)sin(β)  # sun frame theta-direction, not scaled
    a3 = -(C1*cos(α)+C2)*sin(α)cos(β)  # sun frame h-direction, not scaled
    a_SRP = cos(α)*sc.areaParam * G0/d^2 * @SVector [a1; a2; a3]  # aSRP expressed in sun frame
    return a_SRP
end

"""
Function: aSRP_orbitFrame
Description:
    Calculating the SRP acceleration on the solar sail in the sc-centered orbit (LVLH) frame
INPUTS:
spacecraftState: current Keplerian orbital elements of s/c
u: control vector containing:
    α: cone angle (radians)
    β: clock angle (radians)
sc: solar sail struct 
trueAnom_earth: earth current true anomaly (to calculate SRP at)
eph: ephemeride struct
OUTPUTS:
aSRP_oFrame: acceleration vector in the orbit frame [km/s^2]
"""
function aSRP_orbitFrame(spacecraftState, u, sc::basicSolarSail, trueAnom_earth, eph::TwoBodyEphemeride)
    # unpack current state
    method = sc.method
    inc = spacecraftState[3]; 
    argPer = spacecraftState[4];
    trueAnom = spacecraftState[6];
    a_S = aSRP(u, sc, eph, trueAnom_earth)

    # Now, create rotation matrix from sun frame to orbit frame
    # R_S_0 = R3(trueAnom + argPer)*R1(inc)*R3(lambda)
    if method == :Oguri
        lambda = spacecraftState[5]
    elseif method == :Kep
        lambda = spacecraftState[5]-trueAnom_earth
    end
    R_S_O = hill_to_orbit_transform(inc, argPer, lambda, trueAnom)
    aSRP_oFrame = R_S_O * a_S
    return aSRP_oFrame
end


function gauss_variational_eqn!(dx, x, params::QLawParams, t)
    #= 
    Gauss's Variational Equations for the evolution of orbital parameters over time, subject to acceleration
    INPUTS: ALL ANGLES IN STATE VECTOR SHOULD BE IN RADIANS

    OUTPUTS:
    =#
    # unpack parameters:
    mu_sun = 1.327E11
    mu = params.mu
    sc = params.sc
    u = @SVector [params.alpha; params.beta]  # control inputs, alpha and beta
    eph = params.eph  # ephemeride (contains start date)
    epoch = eph.t0 #start date of propagation , et
    method = sc.method

    #Unpack state vector:
    a = x[1]
    e = x[2]
    inc = x[3]
    argPer = x[4]
    lambda = x[5]
    trueAnom = x[6]

    # If eccentricity or inclination would be zero, hold them slightly above (feels crude but necessary)
    if e <= 1.0E-5
        e = 1.0E-5
    end
    if inc <= 1.0E-5
        inc = 1.0E-5
    end
    state = @SVector [a, e, inc, argPer, lambda, trueAnom]
    # get the acceleration:
    # Compute earth coe at current time
    nue = get_heliocentric_position(eph, params.current_time) # note, args ideally are (eph, eph.t0+t), but leave it like this for now (approx.same result)

    a_SRP_O = aSRP_orbitFrame(state, u, sc, nue, eph);  # SRP acceleration resolved into the O (orbit) frame where O{s/c, er_hat, eθ_hat, eh_hat}

    # Calculate some necessary parameters for Gauss's Variational Equations:
    p = a*(1-e^2)  # s/c semi-latus [km]
    r = p/(1+e*cos(trueAnom)) # s/c radial distance from earth [km]
    h = sqrt(mu*p)  # s/c specific angular momentum [km^2/s]

    #Oguri formulation
    nuDot = sqrt(mu/a^3)*(1+e*cos(trueAnom))^2 / (1-e^2)^(3/2) 

    
    ae = eph.semiMajorAxis
    ee = eph.eccentricity
    nuDot_earth = sqrt(mu_sun/ae^3)*(1+ee*cos(nue))^2 / (1-ee^2)^(3/2) 
    
    if method == :Kep
        f0 = [0;0;0;0;0;nuDot]
    else
        f0 = [0;0;0;0; -nuDot_earth; nuDot]
    end
    F = 
    1/h*[
        2*a^2*e*sin(trueAnom) 2*a^2*p/r 0;
        p*sin(trueAnom) (p+r)*cos(trueAnom)+r*e 0;
        0 0 r*cos(trueAnom+argPer);
        -p*cos(trueAnom)/e (p+r)*sin(trueAnom)/e -r*sin(trueAnom+argPer)/tan(inc);
        0 0 r*sin(trueAnom+argPer)/sin(inc);
        p*cos(trueAnom)/e -(p+r)*sin(trueAnom)/e 0;
    ]

    dx[1:6] .= f0 + F*(a_SRP_O); 

    #dx[1:6] .= f0 + F*[0; 1E-7; 0] # uncomment to test EOM
end

"""
Function: QLawEOM
Description:
    Equation of motion with QLaw-calculated control angles
Inputs: 
    dx: derivative vector (for in-place from)
    x: state vector (Keplerian orbital elements)
    params: QLawParams
    t: time [s]
"""
function QLawEOM!(dx, x, params::QLawParams, t)
    # Unpack
    eph = params.eph # ephemeride struct containing earth's helio position at epoch and epoch time in ephemeris time
    mu = params.mu
    a = x[1]
    e = x[2]
    inc = x[3]
    argPer = x[4]
    lambda = x[5]
    trueAnom = x[6]

    #Check for NaN/Infs in time
    if isnan(t) || isinf(t)
        error("Time is $t")
    end

    # Update params current time
    params.current_time = t # current time in ephemeris time [s]

    # Update Earth Position
    nue = get_heliocentric_position(eph, t)

    # Compute Control:
    alphastar, betastar, dQdx = compute_control(x, params)

    # Check for state convergence
    # Stopping criteria
    # Compute targeting error
    aerr                = params.Woe[1]*abs(x[1] - params.oet[1]) - params.oeTols[1]
    eerr                = params.Woe[2]*abs(x[2] - params.oet[2]) - params.oeTols[2]
    ierr                = params.Woe[3]*abs(x[3] - params.oet[3]) - params.oeTols[3]
    ωerr                = params.Woe[4]*abs(acos(cos(x[4] - params.oet[4]))) - params.oeTols[4]
    λerr                = params.Woe[5]*abs(acos(cos(x[5] - params.oet[5]))) - params.oeTols[5]
    targError           = @SVector [aerr, eerr, ierr, ωerr, λerr]
    if maximum(targError) <= 0
            params.terminate = true
    end


    u = @SVector [alphastar; betastar]
    ## Compute derivatives based on control:
    a_SRP_O = aSRP_orbitFrame(x, u, sc, nue, eph);  # SRP acceleration resolved into the O (orbit) frame where O{s/c, er_hat, eθ_hat, eh_hat}
    p = a*(1-e^2)  # s/c semi-latus [km]
    r = p/(1+e*cos(trueAnom)) # s/c radial distance from earth [km]
    h = sqrt(mu*p)  # s/c specific angular momentum [km^2/s]
    nuDot = sqrt(mu/a^3)*(1+e*cos(trueAnom))^2 / (1-e^2)^(3/2) # s/c ballistic evolution
    ae = eph.semiMajorAxis
    ee = eph.eccentricity
    mu_sun = params.mu_sun
    nuDot_earth = sqrt(mu_sun/ae^3)*(1+ee*cos(nue))^2 / (1-ee^2)^(3/2) 
    f0 = [0;0;0;0; -nuDot_earth; nuDot]
    F = 
    1/h*[
        2*a^2*e*sin(trueAnom) 2*a^2*p/r 0;
        p*sin(trueAnom) (p+r)*cos(trueAnom)+r*e 0;
        0 0 r*cos(trueAnom+argPer);
        -p*cos(trueAnom)/e (p+r)*sin(trueAnom)/e -r*sin(trueAnom+argPer)/tan(inc);
        0 0 r*sin(trueAnom+argPer)/sin(inc);
        p*cos(trueAnom)/e -(p+r)*sin(trueAnom)/e 0;
    ]

    
    # Check Lyapunov function [dQdt should be negative]
    dQdt = dQdx'*(f0 + F*(a_SRP_O))
    if  dQdt > 0.0
        @infiltrate false # for debugging
        alphastar = pi/2 # if dQdt>0, zero out SRP acceleration from the sail by setting alpha=90deg
        a_SRP_O = aSRP_orbitFrame(x, [alphastar; betastar], sc, nue, eph);
        #error("Lyapunov condition (dQdt<0) not met at t=$t. dQdt at this time: $dQdt")
    end
    dx[1:6] .= f0 + F*(a_SRP_O); 
    @infiltrate false # this is a good point for debugging, set false->true to turn on breakpoint
end
