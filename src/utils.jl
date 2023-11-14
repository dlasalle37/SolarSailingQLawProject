using LinearAlgebra
# Some utility functions, not specifically associated with anything
function shadow_function()
    #= Formulation can be found sure at https://ai-solutions.com/_help_Files/spacecraft_percentshadow_0_nanosecond.htm\
    Inputs: (edit once defined)
        - Spacecraft position description 
        - some time variable to place the Sun wrt Earth
    Outputs: 
        - k: shading coefficient ∈ [0, 1], 0 for no sunlight, 1 for direct sunlight
    Usage:
        - To be used as a multiplier on the acceleration vector of a solar sailing spacecraft, based on the instantaneous position and time
    =#
    percentShadow = 0; # 1 for fully shaded, 0 for no shade, REPLACE LATER W/ ACTUAL FORMULATION


    k = 1-percentShadow
    return k
end



function rv2coe(r::Vector, v::Vector, mu::Float64) 
    #=
    Calculate the 6 classical orbital elements from r, v, mu (under the two body problem)
    INPUTS:
        r: inertial position vector (km)
        v: inertial velocity vector (km/s)
        mu: gravitational parameter of central body (km^3/s^2)
    OUPUTS:
        a: semi major axis (km)
        e: eccentricity
        inc: inclination (radians)
        RAAN: right-ascension (radians)
        argPer: argument of perigee (radians)
        trueAnom: true anomaly (radians)
    =#
    r_norm = norm(r); #scalar position
    v_norm = norm(v); #scalar velocity

    n_z = [0, 0, 1]; #inertial z unit vector
    n_x = [1, 0, 0]; #inertial x unit vector

    #Semi-Major
    #From vis-viva:
    a = (2/r_norm - v_norm^2/mu)^-1;

    #Eccentricity
    h_vec = cross(r,v); #angular momentum vector
    e_vec = 1/mu * (cross(v,h_vec)) - 1/r_norm * r; #eccentricity vector
    e = norm(e_vec); #eccentricity

    #true anomaly
    nu = acos(dot(r, e_vec) / (r_norm * e)); # true anomaly
    #need to check condition
    #if dot(r,v)>0, 0<nu<180 (no change needed)
    #if dot(r,v)<0, 180<nu<360

    if dot(r,v)<0
        nu = 2*pi - nu;
    end

    #Inclination
    h3  = h_vec[3]; #z component of h
    h = norm(h_vec); #scalar ang. momentum
    i = acos(h3/h); 

    #RAAN
    h = norm(h_vec);
    e_h = h_vec / h; #e_h vector from eq (115) 
    e_n = cross(n_z, e_h);
    en1 = dot(e_n, n_x); #grab first (x) component of e_n
    RAAN = acos(en1/norm(e_n)); 

    #conditional for second component of e_n being negative
    if e_n[2] < 0
        RAAN = 2*pi - RAAN;
    end

    #%argument of perigee
    e_e = e_vec / e; #normalized eccentricity vector in inertial coords
    w = acos(dot(e_n,e_e) / (norm(e_n) * norm(e_e)));
    #conditional:
    if e_e[3] < 0
        w = 2*pi - w;
    end

    coe = (a, e, i, w, RAAN, nu)
    return coe
    
end

function coe2rv(a::Float64, e::Float64, i::Float64, w::Float64, RAAN::Float64, nu::Float64, mu::Float64)
    #=
    Calculate the inertial position and velocity vectors from Keplerian classical orbital elements
    INPUTS:
        a: semi major axis (km)
        e: eccentricity
        inc: inclination (radians)
        RAAN: right-ascension (radians)
        argPer: argument of perigee (radians)
        trueAnom: true anomaly (radians)
        mu: gravitational parameter of central body (km^3/s^2)
    Outputs:
        r: inertial position vector (km)
        v: inertial velocity vector (km/s)
    =#
    p = a*(1-e^2); #semi-latus rectum
    r = p / (1+e*cos(nu)); #scalar r from trajectory

    r_p = [r*cos(nu); r*sin(nu); 0]; #perifocal position vector


    v_p = [-sqrt(mu/p) * sin(nu); sqrt(mu/p)*(e+cos(nu)); 0]; #perifocal v

    #need to transform r_p and v_p  to inertial
    #use p->n transfer matrix, A. (eq126):

    A11 = cos(RAAN)*cos(w) - sin(RAAN)*sin(w)*cos(i);
    A12 = -cos(RAAN)*sin(w) - sin(RAAN)* cos(w)*cos(i);
    A13 = sin(RAAN)*sin(i);
    A21 = sin(RAAN)*cos(w) + cos(RAAN)*sin(w)*cos(i);
    A22 = -sin(RAAN)*sin(w) + cos(RAAN)*cos(w)*cos(i);
    A23 = -cos(RAAN)*sin(i);
    A31 = sin(w)*sin(i);
    A32 = cos(w)*sin(i);
    A33 = cos(i);

    A = [A11 A12 A13;
        A21 A22 A23;
        A31 A32 A33];
    #inertial vectors:

    r = A*r_p; 
    v = A*v_p;
    return r, v
end

function get_sunlight_direction(et)
    #= 
    Calls the SPICE toolkit to calculate the position of the earth WRT the sun expressed in the J2000 Frame.
    Since this vector points to the Earth, this will give the sunlight direction and the radial distance to the earth
    
    Inputs: 
        et: seconds past J2000 date
    Outputs:
        sHat: sunlight direction unit vector in j2k frame
        rSE: radial distance of the sun to the earth (this is approximately equal to sun-spacecraft distance for near-earth satellites)
    =#

    (X0, ~) = spkez(399, et, "J2000", "none", 10)
    rES_j2k = X0[1:3]
    rES = norm(rES_j2k)
    sr = rES_j2k/rES

    return sr, rES

end

function get_heliocentric_position(eph::TwoBodyEphemeride, tf; tol=1.0E-6)
    #=
    INPUTS: 
    eph: ephemeride struct containing the initial time and position info of earth in helio orbit
    tf: current time, seconds past j2000

    OUTPUTS:
        trueAnom: current position of eph.targID in its heliocentric orbit
    =#
    mu = SUN_MU
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

function distance_to_sun(eph::TwoBodyEphemeride, nu::Float64)
    #=
    Calculates the planet-sun distance based on earth's current true anomaly
    INPUTS:
        eph: TwoBodyEphemeride with target as planet and observer as sun
        nu: true anomaly to calculate distance at [radians]
    OUTPUTS:
        dist: distance to sun [km]
    =#
    a = eph.semiMajorAxis
    e = eph.eccentricity

    dist = a*(1-e^2)/(1+e*cos(nu))  # trajectory equation
    return dist
end

"""
hill_to_orbit_transform: A function to construct the rotation matrix from the hill frame to the orbit frame
Notes: 
    - This is a general rotation matrix defined in 'Solar Sailing Q-Law for Planetocentric, Many-Revolution Sail Orbit Transfers' by Oguri & McMahon
Inputs:
    inc: inclination [rad]
    ape: argument of periapsis [rad]
    lam: longitude of ascending node from Hill frame i_r vector [rad]
    tru: true anomaly [rad]
"""
function hill_to_orbit_transform(inc, ape, lam, tru)
    R = @SArray [  # Rotation from hill to orbit frame
            cos(ape + tru)*cos(lam) - sin(ape + tru)*cos(inc)*sin(lam) cos(ape + tru)*sin(lam) + sin(ape + tru)*cos(inc)*cos(lam) sin(ape + tru)*sin(inc);
             -sin(ape + tru)*cos(lam) - cos(ape + tru)*cos(inc)*sin(lam) cos(ape + tru)*cos(inc)*cos(lam)-sin(ape + tru)*sin(lam) cos(ape + tru)*sin(inc);
             sin(inc)*sin(lam) -cos(lam)*sin(inc) cos(inc)
            ]
    return R
end