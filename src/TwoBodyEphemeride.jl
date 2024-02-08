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

"""
twoBodyEarthEphemeride: Constructor to create an ephemeride for Earth's Heliocentric orbit about the sun
Inputs:
    t0: Initial time (in ephemeris time)
    tf: final time (in ephemeris time)
    targID: SPICE-defined id [default 399 for Earth]
    obsID: SPICE-defined id [default 10 for Sun]
"""
function twoBodyEarthEphemeride(t0, tf; targID=399, obsID=10)  #default observer is sun, default target is earth
    (tA_i,~) = getInitialState(t0, targID) # get heliocentric trueanomaly at t0
    ecc = EARTH_HELIOC_ECC
    sma = EARTH_HELIOC_SEMIMAJOR

    TwoBodyEphemeride(t0, tf, targID, obsID, tA_i, sma, ecc)
end

"""
Called in constructor for earth ephemeride
NOTES: Specific to heliocentric orbits. Call spice to find the initial true anomaly of a celestial body about the sun
INPUTS:
    t0: initial epoch in ephemeris time
OUTPUTS:
    trueAnom: true anomaly [radians]
    dist: sun->target body [km]
"""
function getInitialState(t0, targID)
    (X1, ~) = spkez(targID, t0, "ECLIPJ2000", "none", 10);
    r_vec = X1[1:3]; v_vec = X1[4:6]  # split state vector X1 into r, v
    dist = norm(r_vec)
    coe = rv2coe(r_vec, v_vec, SUN_MU)  # keplerian orbital elements in its heliocentric orbit

    trueAnom = coe[6]  # pull out true anomaly [radians]
    return  trueAnom, dist
end

"""
get_heliocentric_position: Solve kepler's equation for ellipses to calculate the position of a planet on its
heliocentric orbit (true anomaly)
INPUTS: 
eph: ephemeride struct containing the initial time and position info of earth in helio orbit
tf: current time in ephemeris time

OUTPUTS:
    trueAnom: current position of eph.targID in its heliocentric orbit [rad]
"""
function get_heliocentric_position(eph::TwoBodyEphemeride, tf; tol=1.0E-6)
    mu = SUN_MU
    t0 = eph.t0;
    t = tf-t0 # difference between current and initial ephemeris times [seconds]
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
    numLoops = 0
    while !done
        nextGuess = Ek + (M0 + M - Ek + e*sin(Ek))/(1-e*cos(Ek))
        err = abs(nextGuess-Ek)

        if err <= tol
            Ek = nextGuess
            done = true
        else
            Ek = nextGuess
            done = false
            if numLoops >= 100
                
                done = true
                println(eph)
                println("tf: $tf")
                error("Earth position calculation nonconvergence")
                
            end
        end
        numLoops+=1
    end
    E=Ek
    νf = 2*atan(sqrt((1+e)/(1-e)) * tan(E/2))
    if νf < 0 # putting νf on range 0-2pi
        νf += 2*pi
    end

    return νf
end

"""
Calculates the planet-sun distance based on earth's current true anomaly.
Just applies the trajectory equation, based on nominal values of semi-major(a) and eccentricity(e).

Later, it may be possible to implement oscillating values for a and e, just need to turn the TwoBodyEphemeride struct to mutable, and continually update
a and e values.
INPUTS:
    eph: TwoBodyEphemeride with target as planet and observer as sun
    nu: true anomaly to calculate distance at [radians]
OUTPUTS:
    dist: distance to sun [km]
"""
function distance_to_sun(eph::TwoBodyEphemeride, nu::Float64)
    a = eph.semiMajorAxis
    e = eph.eccentricity

    dist = a*(1-e^2)/(1+e*cos(nu))  # trajectory equation
    return dist
end