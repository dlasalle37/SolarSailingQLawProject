
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


"""
Calculate the 6 classical keplerian orbital elements from r, v, mu (under the two body problem)
INPUTS:
    r: inertial position vector (km)
    v: inertial velocity vector (km/s)
    mu: gravitational parameter of central body (km^3/s^2)
OUPUTS:
    coe: tuple of 6 Keplerian elements, which are:
        a: semi major axis (km)
        e: eccentricity
        inc: inclination (radians)
        RAAN: right-ascension (radians)
        argPer: argument of perigee (radians)
        trueAnom: true anomaly (radians)
"""
function rv2coe(r, v, mu) 
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

"""
coe2rv: Calculate the inertial position and velocity vectors from Keplerian classical orbital elements
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
"""
function coe2rv(a, e, i, w, RAAN, nu, mu)
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

"""
Calls the SPICE toolkit to calculate the position of the earth WRT the sun expressed in the J2000 Frame.
Since this vector points to the Earth, this will give the sunlight direction and the radial distance to the earth

Inputs: 
    et: seconds past J2000 date
Outputs:
    sHat: sunlight direction unit vector in j2k frame
    rES: radial distance of the sun to the earth (this is approximately equal to sun-spacecraft distance for near-earth satellites)
"""
function get_sunlight_direction(et)
    (X0, ~) = spkez(399, et, "J2000", "none", 10)
    rES_j2k = X0[1:3]
    rES = norm(rES_j2k)
    sr = rES_j2k/rES # sunlight direction vector expressed in j2000 frame

    return sr, rES

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


"""
function get_gm(id): get the gravitational parameter (GM, or mu) of a celestial body by passing in its SPICE ID
    INPUTS: id: SPICE ID of body
    OUTPUTS: mu: GM of body given by [id]
    
    dependencies: Info is pulled from GRAV_PARAMS dictionary in constants.jl
"""
function get_gm(id)
    mu = GRAV_PARAMS[id]
    return mu
end