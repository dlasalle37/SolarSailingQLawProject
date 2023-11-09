using LinearAlgebra
using SPICE

struct basicSolarSail
    areaParam::Float64
    ρ::Float64
    s::Float64
    ϵf::Float64
    ϵb::Float64
    Bf::Float64
    Bb::Float64
    C::Vector{Float64} # contains all 3 (C1, C2, C3) condensed parameters
    accelerationModel::Symbol  # acceleration model to be used; default is FPOF for flat-plate optical force (Model used in Oguri)
end

# constructor
function basicSolarSail(; areaParam=0.0025/100, ρ=1, s=1, ϵf=0.1, ϵb=0.1, Bf=2/3, Bb=2/3, accelerationModel=:FPOF)

    #Checking bounds
    if !(0<=ρ<=1) || !(0<=s<=1)
        error("Sail parameters ρ and/or s out of bounds [0, 1]")
    end
    (C1, C2, C3) = calculate_sail_constants(ρ, s, ϵf, ϵb, Bf, Bb)
    condensedParams = [C1,C2,C3]

    # Create struct
    basicSolarSail(areaParam, ρ, s, ϵf, ϵb, Bf, Bb, condensedParams, accelerationModel)
end

function calculate_sail_constants(ρ, s, ϵf, ϵb, Bf, Bb)
    #= 
    This function executes in the constructor.
    Uses the sail parameters to calculate the condensed parameters C1, C2, C3
    Formulation discussed in Solar Sailing Primer Vector Theory: Indirect Trajectory
    Optimization with Practical Mission Considerations by Oguri, McMahon, Lantoine.
    =#

    C1 = 2*ρ*s
    C3 = 1-ρ*s
    C2 = Bf*ρ*(1-s)+(1-ρ)*(Bf*ϵf-Bb*ϵb)/(ϵf+ϵb)

    if (C1==0) && (C2==0)
        error("Check sail parameters, C1 and C2 cannot both be zero")
    elseif !(0<=C1<=2) || !(0<=C2<=1)
        error("Values out of bounds. C1∈[0,2] and C2∈[0,1]")
    end
    return C1, C2, C3
end

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
    G0 = 1.02E14
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

function aSRP_orbitFrame(spacecraftState, α::Float64, β::Float64, sc::basicSolarSail, time::Float64, eph::TwoBodyEphemeride; method=:Oguri)
    #=
    INPUTS:
        spacecraftState: current Keplerian orbital elements of s/c
        α: cone angle (radians)
        β: clock angle (radians)
        sc: solar sail struct 
        time: time past j2000 to calculate SRP
        mu: grav parameter of central body [km^3/s^2]
        method: :Oguri (default) for a state vector of [a e inc argPer lambda, trueAnom]
                :Kep for state vector of [a e inc argPer RAAN, trueAnom]
    OUTPUTS:
        aSRP_oFrame: acceleration vector in the orbit frame [km/s^2]
        earth_coe: keplerian orbital elements of earth at current time (saving space)
    =#
    # unpack current state
    a = spacecraftState[1]; 
    e = spacecraftState[2]; 
    inc = spacecraftState[3]; 
    argPer = spacecraftState[4];
    trueAnom = spacecraftState[6];

    G0 = 1.02E14 # solar flux constant at earth
    trueAnom_earth = get_heliocentric_position(eph, time)
    d = distance_to_sun(eph, trueAnom_earth)
    C1 = sc.C[1]; C2 = sc.C[2]; C3 = sc.C[3]
    a1 = C1*(cos(α))^2+C2*cos(α)+C3  # sun frame r-direction, not scaled
    a2 = -(C1*cos(α)+C2)*sin(α)sin(β)  # sun frame theta-direction, not scaled
    a3 = -(C1*cos(α)+C2)*sin(α)cos(β)  # sun frame h-direction, not scaled
    a_S = cos(α)*sc.areaParam * G0/d^2 * @SVector [a1; a2; a3]  # aSRP expressed in sun frame

    # Now, create rotation matrix from sun frame to orbit frame
    # R_S_0 = R3(trueAnom + argPer)*R1(inc)*R3(lambda)
    if method == :Oguri
        lambda = spacecraftState[5]
    elseif method == :Kep
        lambda = spacecraftState[5]-trueAnom_earth
    end
    ang = trueAnom + argPer;
    R3a = @SMatrix [cos(ang) sin(ang) 0; -sin(ang) cos(ang) 0; 0 0 1]
    R1 = @SMatrix [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)]
    R3b = @SMatrix [cos(lambda) sin(lambda) 0; -sin(lambda) cos(lambda) 0; 0 0 1]
    R_S_O = R3a*R1*R3b
    aSRP_oFrame = R_S_O * a_S
    return aSRP_oFrame
end


function gauss_variational_eqn!(dx, x, params, t)
        #= 
    Gauss's Variational Equations for the evolution of orbital parameters over time, subject to acceleration
    INPUTS: ALL ANGLES IN STATE VECTOR SHOULD BE IN RADIANS

    OUTPUTS:
    =#
    # unpack parameters:
    mu_sun = 1.327E11
    mu = params[1]
    sc = params[2]
    u = params[3]  # control inputs, alpha and beta
    eph = params[4]  # ephemeride (contains start date)
    epoch = eph.t0 #start date of propagation , et
    α = u[1]  # cone angle
    β = u[2]  # clock angle

    #Unpack state vector:
    a = x[1]
    e = x[2]
    inc = x[3]
    argPer = x[4]
    lambda = x[5]
    trueAnom = x[6]

    # If eccentricity or inclination would be zero, hold them slightly above (feels crude but necessary)
    if e <= 1.0E-4
        e = 1.0E-4
    end
    if inc <= 1.0E-4
        inc = 1.0E-4
    end
    state = (a, e, inc, argPer, lambda, trueAnom)
    # get the acceleration:
    a_SRP_O = aSRP_orbitFrame(state, α, β, sc, epoch+t, eph);  # SRP acceleration resolved into the O (orbit) frame where O{s/c, er_hat, eθ_hat, eh_hat}

    # Calculate some necessary parameters for Gauss's Variational Equations:
    p = a*(1-e^2)  # s/c semi-latus [km]
    r = p/(1+e*cos(trueAnom)) # s/c radial distance from earth [km]
    h = sqrt(mu*p)  # s/c specific angular momentum [km^2/s]

    #Oguri formulation
    nuDot = sqrt(mu/a^3)*(1+e*cos(trueAnom))^2 / (1-e^2)^(3/2) 

    # Compute earth coe at current time
    nue = get_heliocentric_position(eph, epoch+t)
    ae = eph.semiMajorAxis
    ee = eph.eccentricity
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

    dx[1:6] .= f0 + F*(a_SRP_O); 

    #dx[1:6] .= f0 + F*[0; 1E-7; 0] # uncomment to test EOM
end

function gauss_variational_eqn_keplerian!(dx, x, params, t)
    #= 
    Gauss's Variational Equations for the evolution of orbital parameters over time, subject to acceleration
    INPUTS: ALL ANGLES IN STATE VECTOR SHOULD BE IN RADIANS

    OUTPUTS:
    =#
    # unpack parameters:
    mu_sun = 1.327E11
    mu = params[1]
    sc = params[2]
    u = params[3]  # control inputs, alpha and beta
    eph = params[4]  # ephemeride (contains start date)
    epoch = eph.t0 #start date of propagation , et
    α = u[1]  # cone angle
    β = u[2]  # clock angle

    #Unpack state vector:
    a = x[1]
    e = x[2]
    inc = x[3]
    argPer = x[4]
    RAAN = x[5]
    trueAnom = x[6]

    # If eccentricity or inclination would be zero, hold them slightly above (feels crude but necessary)
    if e <= 1.0E-4
        e = 1.0E-4
    end
    if inc <= 1.0E-4
        inc = 1.0E-4
    end
    state = (a, e, inc, argPer, RAAN, trueAnom)
    # get the acceleration:
    a_SRP_O = aSRP_orbitFrame(state, α, β, sc, epoch+t, eph, method=:Kep);  # SRP acceleration resolved into the O (orbit) frame where O{s/c, er_hat, eθ_hat, eh_hat}

    # Calculate some necessary parameters for Gauss's Variational Equations:
    p = a*(1-e^2)  # s/c semi-latus [km]
    r = p/(1+e*cos(trueAnom)) # s/c radial distance from earth [km]
    h = sqrt(mu*p)  # s/c specific angular momentum [km^2/s]

    Akep = [0;0;0;0;0; h/r^2]

    Bkep = 1/h*[
        2*a^2*e*sin(trueAnom) 2*a^2*p/r 0;
        p*sin(trueAnom) (p+r)*cos(trueAnom)+r*e 0;
        0 0 r*cos(trueAnom+argPer);
        -p*cos(trueAnom)/e (p+r)*sin(trueAnom)/e -r*sin(trueAnom+argPer)/tan(inc);
        0 0 r*sin(trueAnom+argPer)/sin(inc);
        p*cos(trueAnom)/e -(p+r)*sin(trueAnom)/e 0
        ]
    
    dx[1:6] = Akep + Bkep*a_SRP_O
end    