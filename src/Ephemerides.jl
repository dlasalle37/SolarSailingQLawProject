# Ephemeride struct
struct TwoBodyEphemeride
    t0::Float64
    tf::Float64
    targID::Int64
    obsID::Int64

    # Initial orbital parameters
    trueAnom_initial::Float64
    semiMajorAxis::Float64
    eccentricity::Float64
end

#constructor
function twoBodyEarthEphemeride(t0, tf; targID=399, obsID=10)  #default observer is sun, default target is earth
    (tA_i,~) = getInitialState(t0, targID)
    ecc = EARTH_HELIOC_ECC
    sma = EARTH_HELIOC_SEMIMAJOR

    TwoBodyEphemeride(t0, tf, targID, obsID, tA_i, sma, ecc)
end

function getInitialState(t0, targID)
    #= Called in constructor for earth ephemeride
    NOTES: Specific to heliocentric orbits
    INPUTS:
        t0: initial epoch in time past j2000 [s]
    OUTPUTS:
        trueAnom: true anomaly [radians]
        dist: sun->target body
    =#

    (X1, ~) = spkez(targID, t0, "ECLIPJ2000", "none", 10);
    r_vec = X1[1:3]; v_vec = X1[4:6]  # split state vector X1 into r, v
    dist = norm(r_vec)
    coe = rv2coe(r_vec, v_vec, 1.327E11)  # keplerian orbital elements in its heliocentric orbit

    trueAnom = coe[6]  # pull out true anomaly [radians]
    return  trueAnom, dist
end

function earth_heliocentric_position(eph::TwoBodyEphemeride, tf; tol=1.0E-6)
    #=
    INPUTS: 
    eph: ephemeride struct containing the initial time and position info of earth in helio orbit
    tf: current time, seconds past j2000

    OUTPUTS:
        trueAnom: current position of earth in its heliocentric orbit
    =#
    mu = MU_SUN
    t0 = eph.t0;
    t = tf-t0
    νi = eph.trueAnom_initial
    a = eph.semiMajorAxis
    e = eph.eccentricity

    # use keplers equation for ellipses to calculate current position
    n = sqrt(mu/a^3)
    M = n*t  # current mean anomaly to calculate to
    E0 = 2*atan(sqrt((1-e)/(1+e))*tan(νi/2))  # initial eccentric anomaly
    M0 = E0 - e*sin(E0)

    #Use newton's method to calculate E at t;
    done = false
    Ek = M  # initial guess
    while !done
        nextGuess = Ek + (M0 + M - Ek + e*sin(Ek))/(1-e*cos(Ek))
        err = abs(nextGuess-Ek)

        if err <= tol
            Ek = nextGuess
            done = true
        else
            Ek = nextGuess
            done = false
        end
    end
    E=Ek
    νf = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2))
    return νf
end