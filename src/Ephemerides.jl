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

