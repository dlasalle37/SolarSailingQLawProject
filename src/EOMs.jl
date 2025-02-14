"""
    keplerian_varational_equations(a, e, i, ω, Ω, θ, mu)
Stacked Gaussian variational equations for the orbital element set [a, e, i, ω, Ω, θ] (defined below)

# Inputs:
    a: semi-major axis [km]
    e: eccentricity
    i: inclination
    ω: argument of perigee [rad]
    Ω: right ascencion of ascending node [rad]
    θ: true anomaly [rad]
# Outputs:
    f0: 6x1 ballistic evolution vector 
    F: 6x3 Variational matrix to be multiplied with an acceleration (3x1) in the sc-centered orbit frame
"""
function keplerian_variational_equations(a, e, i, ω, Ω, θ, mu)
    p = a*(1-e^2)  # s/c semi-latus [km]
    r = p/(1+e*cos(θ)) # s/c radial distance from earth [km]
    h = sqrt(mu*p)  # s/c specific angular momentum [km^2/s]
    θ_dot = sqrt(mu/a^3)*(1+e*cos(θ))^2 / (1-e^2)^(3/2)
    f0 = @SVector [0;0;0;0;0; θ_dot]
    F = 
    @SArray [
        2*a^2*e*sin(θ) 2*a^2*p/r 0;
        p*sin(θ) (p+r)*cos(θ)+r*e 0;
        0 0 r*cos(θ+ω);
        -p*cos(θ)/e (p+r)*sin(θ)/e -r*sin(θ+ω)/tan(i);
        0 0 r*sin(θ+ω)/sin(i);
        p*cos(θ)/e -(p+r)*sin(θ)/e 0;
    ]
    F = F*1/h 

    return f0, F
end

"""
    augmented_keplerian_varaitions:
Stacked Gaussian variational equations for the orbital element set [a, e, i, ω, λ, θ] (defined below)

# INPUTS:
    a: semi-major axis [km]
    e: eccentricity
    i: inclination
    ω: argument of perigee [rad]
    λ: ascencion of ascending node measured from the ir_hat vector (vector from sun to central body(eg. Earth) direction) [rad]
    θ: true anomaly [rad]
    θ_dot_body: ballistic evolution term of central body (such as earth) on its heliocentric orbit [rad/s]
    mu: gravitational parameter of central body [km^3/s^2]

# OUTPUTS:
    f0: 6x1 ballistic evolution vector 
    F: 6x3 Variational matrix to be multiplied with an acceleration (3x1) in the sc-centered orbit frame
"""
function augmented_keplerian_varaitions(a, e, i, ω, λ, θ, θ_dot_body, mu)


    p = a*(1-e^2)  # s/c semi-latus [km]
    r = p/(1+e*cos(θ)) # s/c radial distance from earth [km]
    h = sqrt(mu*p)  # s/c specific angular momentum [km^2/s]
    θ_dot = sqrt(mu/a^3)*(1+e*cos(θ))^2 / (1-e^2)^(3/2)
    f0 = @SVector [0;0;0;0; -θ_dot_body; θ_dot]
    F = 
    @SArray [
        2*a^2*e*sin(θ) 2*a^2*p/r 0;
        p*sin(θ) (p+r)*cos(θ)+r*e 0;
        0 0 r*cos(θ+ω);
        -p*cos(θ)/e (p+r)*sin(θ)/e -r*sin(θ+ω)/tan(i);
        0 0 r*sin(θ+ω)/sin(i);
        p*cos(θ)/e -(p+r)*sin(θ)/e 0;
    ]
    F = F*1/h 

    return f0, F
end

"""
    aSRP
calculating SRP acceleration in the Hill (Sun centered) frame

# Notes:
    - SRP acceleration vector is calculated in the sun-centered hill frame
# inputs:
    u: control variables [alpha, beta]
    sc: basicSolarSail
    eph: ephemeride
    t: ephemeris time [s]

# outputs:
    aSRP: 3x1 SRP acceleration vector in Hill (sun-centered) frame
"""
function aSRP(u, sc::basicSolarSail, eph::Ephemeride, t)
    # Unpack Inputs:
    α = u[1]
    β = u[2]
    C1 = sc.C[1]
    C2 = sc.C[2]
    C3 = sc.C[3]
    C1 = sc.C[1]; C2 = sc.C[2]; C3 = sc.C[3]

    # Ephemeride-based calculations
    ephState = getState(eph, t)
    d = norm(view(ephState, 1:3)) # distance

    G0 = get_solar_flux(eph.targ) # solar flux constant at earth (or other target of ephemeride)[kgkm/s^2]
    a1 = C1*cos(α)^2+C2*cos(α)+C3  # sun frame r-direction, not scaled
    a2 = -(C1*cos(α)+C2)*sin(α)sin(β)  # sun frame theta-direction, not scaled
    a3 = -(C1*cos(α)+C2)*sin(α)cos(β)  # sun frame h-direction, not scaled
    a_SRP = cos(α)*sc.areaParam * G0/d^2 * @SVector [a1; a2; a3]  # aSRP expressed in sun frame
    #=
    a = α; b = β
    a_SRP = sc.areaParam*G0/d^2 * cos(a)*[
        C1*cos(a)^2+C2*cos(a) + C3;
        -((C1*cos(a)+C2)*sin(a)*sin(b));
        -((C1*cos(a)+C2)*sin(a)*cos(b))
    ]=#

    return a_SRP
end

"""
    aSRP_orbitFrame
# Description:
    Calculating the SRP acceleration on the solar sail in the sc-centered orbit (LVLH) frame
# INPUTS:
spacecraftState: current Keplerian orbital elements of s/c
u: control vector containing:
    α: cone angle (radians)
    β: clock angle (radians)
ps::QLawParams{Oguri}: QLaw parameters containing type info and solar sail struct 
eph: Ephemeride struct (sun should be center)
t: current time in et
# OUTPUTS:
aSRP_oFrame: acceleration vector in the orbit frame [km/s^2]
"""
function aSRP_orbitFrame(spacecraftState, u, ps::QLawParams{Oguri}, eph::Ephemeride, t)
    sc = ps.sc
    # unpack current state
    inc = spacecraftState[3]; 
    argPer = spacecraftState[4];
    trueAnom = spacecraftState[6];
    a_S = aSRP(u, sc, eph, t) # hill frame SRP

    # Now, create rotation matrix from sun frame to orbit frame
    # R_S_0 = R3(trueAnom + argPer)*R1(inc)*R3(lambda)
    lambda = spacecraftState[5]
    R_S_O = hill_to_orbit_transform(inc, argPer, lambda, trueAnom)
    aSRP_oFrame = R_S_O * a_S
    return aSRP_oFrame
end

"""
    aSRP_orbitFrame
# Description:
    Calculating the SRP acceleration on the solar sail in the sc-centered orbit (LVLH) frame
# INPUTS:
spacecraftState: current Keplerian orbital elements of s/c
u: control vector containing:
    α: cone angle (radians)
    β: clock angle (radians)
sc::basicSolarSail{Keplerian}: solar sail struct 
eph: Ephemeride struct (sun should be center)
t: current time in et
# OUTPUTS:
aSRP_oFrame: acceleration vector in the orbit frame [km/s^2]
"""
function aSRP_orbitFrame(spacecraftState, u, ps::QLawParams{Keplerian}, eph::Ephemeride, t)
    sc = ps.sc
    # unpack current state
    inc = spacecraftState[3]; 
    argPer = spacecraftState[4];
    trueAnom = spacecraftState[6];
    a_S = aSRP(u, sc, eph, t) # hill frame SRP

    # Now, create rotation matrix from sun frame to orbit frame
    # R_S_0 = R3(trueAnom + argPer)*R1(inc)*R3(lambda)
    coee = getCOE(eph, t)
    trueAnom_earth = coee[6]
    lambda = spacecraftState[5]-trueAnom_earth

    R_S_O = hill_to_orbit_transform(inc, argPer, lambda, trueAnom)
    aSRP_oFrame = R_S_O * a_S
    return aSRP_oFrame
end




"""
    QLawEOM
# Description:
    Equation of motion with QLaw-calculated control angles
# Inputs: 
    dx: derivative vector (for in-place from)
    x: state vector [semi-major axis, eccentricity, inclination, arg. per., lambda(RAAN Analog from Oguri paper), true anomaly]
    params: QLawParams
    t: time [s]
"""
function QLawEOM!(dx, x, params::QLawParams{Oguri}, t)
    
    # update Time
    # params.current_time = t # current time in ephemeris time [s]

    # Unpack
    eph = params.eph # ephemeride struct containing earth's helio position at epoch and epoch time in ephemeris time
    mu = params.mu
    a = x[1]
    e = x[2]
    inc = x[3]
    ape = x[4]
    lam = x[5]
    tru = x[6]

    #Check for NaN/Infs in time
    if isnan(t) || isinf(t)
        error("Time is $t")
    end

    # Update Earth Position
    coee = getCOE(eph, t)
    nue = coee[6]

    ae = coee[1]
    ee = coee[2]
    mu_sun = params.mu_sun
    nuDot_earth = sqrt(mu_sun/ae^3)*(1+ee*cos(nue))^2 / (1-ee^2)^(3/2) 

    # Get Dynamics at current time
    f0, F = augmented_keplerian_varaitions(a, e, inc, ape, lam, tru, nuDot_earth, mu)

    # Eclipse Check
    RAAN = lam+nue # get the actual keplerian element RAAN from definition of lambda
    (r, v) = coe2rv(a, e, inc, ape, RAAN, tru, mu)
    if !isEclipsed(eph, r, t)
        # Compute Control:
        alphastar, betastar, dQdx = compute_control(t, x, params)
        u = @SVector [alphastar; betastar]

        ## Compute derivatives based on control:
        a_SRP_O = aSRP_orbitFrame(x, u, params, eph, t);  # SRP acceleration resolved into the O (orbit) frame where O{s/c, er_hat, eθ_hat, eh_hat}
        
        # Check Lyapunov function [dQdt should be negative]
        dQdt = dQdx'*(f0 + F*(a_SRP_O))
        if  dQdt > 0.0 # If dQdt> 0, it is implied that there is no control to decrease error
            alphastar = pi/2 # if dQdt>0, zero out SRP acceleration from the sail by setting alpha=90deg
            a_SRP_O = aSRP_orbitFrame(x, [alphastar; betastar], params, eph, t);
            #error("Lyapunov condition (dQdt<0) not met at t=$t. dQdt at this time: $dQdt")
        end

    else
        a_SRP_O = @SVector(zeros(3)) # If isEclipsed is true, SRP is zero
    end

    dx[1:6] .= f0 + F*(a_SRP_O); 

end

"""
    QLawEOMKeplerian
# Description:
    Equation of motion with QLaw-calculated control angles using a traditional keplerian set of orbital elements.
# Inputs: 
    dx: derivative vector (for in-place from)
    x: state vector (Keplerian orbital elements)
    params: QLawParams
    t: time [s]
"""
function QLawEOM!(dx, x, params::QLawParams{Keplerian}, t)
    
    # update Time
    # params.current_time = t # current time in ephemeris time [s]
    # Unpack
    eph = params.eph # ephemeride struct containing earth's helio position at epoch and epoch time in ephemeris time
    mu = params.mu
    a = x[1]
    e = x[2]
    inc = x[3]
    ape = x[4]
    ran = x[5]
    tru = x[6]
    
    #Check for NaN/Infs in time
    if isnan(t) || isinf(t)
        error("Time is $t")
    end

    # Get Dynamics at current time
    f0, F = keplerian_variational_equations(a, e, inc, ape, ran, tru, mu)

    # Eclipse Check
    (r, v) = coe2rv(a, e, inc, ape, ran, tru, mu)
    if !isEclipsed(eph, r, t)
        # Compute Control:
        alphastar, betastar, dQdx = compute_control(t, x, params)
        u = @SVector [alphastar; betastar]

        ## Compute derivatives based on control:
        a_SRP_O = aSRP_orbitFrame(x, u, params, eph, t);  # SRP acceleration resolved into the O (orbit) frame where O{s/c, er_hat, eθ_hat, eh_hat}
        
        # Check Lyapunov function [dQdt should be negative]
        dQdt = dQdx'*(f0 + F*(a_SRP_O))
        if  dQdt > 0.0 # If dQdt> 0, it is implied that there is no control to decrease error
            alphastar = pi/2 # if dQdt>0, zero out SRP acceleration from the sail by setting alpha=90deg
            a_SRP_O = aSRP_orbitFrame(x, [alphastar; betastar], params, eph, t);
            #error("Lyapunov condition (dQdt<0) not met at t=$t. dQdt at this time: $dQdt")
        end

    else
        a_SRP_O = @SVector(zeros(3)) # If isEclipsed is true, SRP is zero
    end
    
    # add in perturbations
    mdl = params.gravity_model
    fs = params.fs
    itrf2gcrf = rotation3(fs, 6, 23, t)[1]; # from ICRF(SPICE ID 6) to GCRF(SPICE ID 23)
    pos_fixed = transpose(itrf2gcrf) * r
    a_perturb_fixed = getFirstPartial(mdl, pos_fixed, false) # get perturbation from gravity in fixed frame
    a_perturb_eci = itrf2gcrf * a_perturb_fixed
    a_perturb_orbit = eci2ric(r,v)*a_perturb_eci
    
    dx[1:6] .= f0 + F*(a_SRP_O+a_perturb_orbit); 
    # dx[1:6] .= f0 + F*(a_SRP_O); 

end

"""
callback_function_error_check: function to be called when creating the callback condition

# INPUTS:
    u: state passed in by condition function
    t: time passed in by condition functiom
    params: QLaw Params containing the weights and tolerances
# OUTPUTS:
    returns true to terminate integration
    returns false to continue integration
"""
function callback_function_error_check(x, t, params::QLawParams)
    aerr                = params.Woe[1]*abs(x[1] - params.oet[1]) - params.oeTols[1]
    eerr                = params.Woe[2]*abs(x[2] - params.oet[2]) - params.oeTols[2]
    ierr                = params.Woe[3]*abs(x[3] - params.oet[3]) - params.oeTols[3]
    ωerr                = params.Woe[4]*abs(acos(cos(x[4] - params.oet[4]))) - params.oeTols[4]
    λerr                = params.Woe[5]*abs(acos(cos(x[5] - params.oet[5]))) - params.oeTols[5]
    targError           = @SVector [aerr, eerr, ierr, ωerr, λerr]
    if maximum(targError) <= 0
        ret = true
    else 
        ret = false
    end
    if x[2] > 1 || x[2] < 0
        ret = true # terminate if eccentricity is out of bounds
    end
    return ret
end

"""
continuous_callback_errorcheck: function to be called when creating a continuous callback

# INPUTS:
    u: state passed in by condition function
    t: time passed in by condition functiom
    params: QLaw Params containing the weights and tolerances
# OUTPUTS:
    ret: maximum of errors, when this reaches zero, the solution has converged
"""
function continuous_callback_errorcheck(x, t, params::QLawParams)
    aerr                = params.Woe[1]*abs(x[1] - params.oet[1]) - params.oeTols[1]
    eerr                = params.Woe[2]*abs(x[2] - params.oet[2]) - params.oeTols[2]
    ierr                = params.Woe[3]*abs(x[3] - params.oet[3]) - params.oeTols[3]
    ωerr                = params.Woe[4]*abs(acos(cos(x[4] - params.oet[4]))) - params.oeTols[4]
    λerr                = params.Woe[5]*abs(acos(cos(x[5] - params.oet[5]))) - params.oeTols[5]
    targError           = @SVector [aerr, eerr, ierr, ωerr, λerr]

    ret = maximum(targError)

    return ret

end
"""
    two_body_eom!(dx, x, mu, t)
2Body equations of motion

# Notes: 
    Only used for plotting initial/target orbits

# INPUTS:
    dx: for inplace-form of diffeq
    x: anonymous state vector:
    mu: central body grav. parameter
    t: anonymous time variable
"""
function two_body_eom!(dx, x, mu, t)
    rvec = x[1:3]
    vvec = x[4:6]

    r = norm(rvec)

    dx[1:6] .= [vvec; -mu/r^3 * rvec]
end

"""
    two_body_eom_perturbed!: 
2Body equations of motion with pertubations due to nonspherical gravity
    
# Notes: 
    Only used for plotting initial/target orbits or testing
    
# INPUTS:
    dx: for inplace-form of diffeq
    x: anonymous state vector:
    ps::Tuple: list of all parameters
    mu: Central body grav. parameter [km/s^2]
    fs::FrameSystem containing GCRF and ITRF frames, as well as a Spacecraft point that can be updated
    mdl::NormalizedGravityModel: Struct containing all defining terms of gravity model
    t: anonymous time variable
"""
function two_body_eom_perturbed!(dx, x, ps::Tuple, t)
    # Set some basic terms
    rvec = x[1:3]
    vvec = x[4:6]

    r = norm(rvec)

    # Load from parameters
    mu = ps[1]
    mdl = ps[2]
    fs = ps[3]


    itrf2gcrf = rotation3(fs, 6, 23, t)[1]  # itrf (id: 6) to gcrf(id:23) rotation matrix (assuming negligible difference between J2000 and GCRF)
    pos_fixed = transpose(itrf2gcrf) * view(x, 1:3)
    a_perturb_fixed = getFirstPartial(mdl, pos_fixed, false) # get perturbation from gravity in fixed frame
    a_perturb = itrf2gcrf * a_perturb_fixed

    a_sum = (-mu/r^3 * rvec) + a_perturb
    dx[1:6] .= [vvec; a_sum]
end

"""
    two_body_eom_perturbed_tb!: 
2Body equations of motion with pertubations due to third body perturbations
    
# Notes: 
    Only used for plotting initial/target orbits or testing
    
# INPUTS:
    dx: for inplace-form of diffeq
    x: anonymous state vector:
    ps::Tuple: list of all parameters
    mu: Central body grav. parameter [km/s^2]
    fs::FrameSystem containing GCRF and ITRF frames, as well as a Spacecraft point that can be updated
    mdl::NormalizedGravityModel: Struct containing all defining terms of gravity model
    t: anonymous time variable
"""
function two_body_eom_perturbed_tb!(dx, x, ps::Tuple, t)
    # Set some basic terms
    rvec = x[1:3]
    vvec = x[4:6]

    r = norm(rvec)

    # Load from parameters
    mu = ps[1]
    eph = ps[2]
    thirdbody = ps[3]

    a_perturb = third_body_perturb(rvec, t, eph, thirdbody)

    a_sum = (-mu/r^3 * rvec) + a_perturb

    dx[1:6] .= [vvec; a_sum]
end

