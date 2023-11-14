"""
Computing Q and dQdx from the state variables and parameters
dQdx is computed numerically with finite difference

Notes:
When calling this function in a script, create a global QLawParams struct first, since I have not yet found a way to include parameters
    into this function (as an argument), since it is used to calculate nuerical derivatives with FiniteDiff
"""
function calculate_Q_derivative(x)

    Q = calculate_Q(x)

    # Calculate derivatives

    return Q
end

"""
The below function is where the actual calculation of Q is done
"""
function calculate_Q(x)

        # Unpack inputs
        a = x[1]
        e = x[2]
        inc = x[3]
        ape = x[4]
        lam = x[5]
        tru = x[6]

        ## Unpack parameters
        # Basic params
        mu = params.mu
        mu_sun = params.mu_sun
        rpmin = params.rp_min
        aesc = params.a_esc
        t = params.current_time

        # Scaling terms
        mpet = params.m_petro
        rpet = params.r_petro
        npet = params.n_petro

        # Targets
        a_t = params.a_t
        e_t = params.e_t
        inc_t = params.inc_t
        ape_t = params.ape_t
        lam_t = params.lam_t

        # Weights
        Wa = params.Wa
        We = params.We
        Winc = params.Winc
        Wape = params.Wape
        Wlam = params.Wlam

        # Ephemeride
        eph = params.eph
    
        ## Create Q
        # Some initial terms
        p = a*(1-e^2) # semi-latus rectum
        h = sqrt(mu*p) # spec. angular momentum

        # Element Selection vectors
        eps_a = [1; 0; 0; 0; 0; 0]
        eps_e = [0; 1; 0; 0; 0; 0] 
        eps_inc = [0; 0; 1; 0; 0; 0]
        eps_ape = [0; 0; 0; 1; 0; 0]
        eps_lam = [0; 0; 0; 0; 1; 0]

        # Penalty Functions
        Wp = params.Wp
        impact_constraint = 1 - a*(1-e)/rpmin
        escape_constraint = a/aesc - 1
        P1 = params.Aimp*exp(params.kimp*impact_constraint)
        P2 = params.Aesc*exp(params.kesc*escape_constraint)
        P = Wp*(P1 + P2)

        # Scaling (a only)
        Sa = (1+((a-a_t)/(mpet*a_t))^npet)^(1/rpet);

        # Distance Functions
        dista = a - a_t;
        diste = e - e_t;
        distinc = inc - inc_t;
        distape = acos(cos(ape - ape_t));
        distlam = acos(cos(lam-lam_t));

        # Sign Functions
        sig_a = sign(dista);
        sig_e = sign(diste);
        sig_inc = sign(distinc);
        sig_ape = sign(distape);
        sig_lam = sign(distlam);

        # Sunlight direction in perifocal frame
        term1 = @SArray     [cos(ape) sin(ape) 0; -sin(ape) cos(ape) 0; 0 0 1];
        term2 = @SArray     [1 0 0; 0 cos(inc) sin(inc); 0 -sin(inc) cos(inc)];
        term3 = @SArray     [cos(lam) sin(lam) 0; -sin(lam) cos(lam) 0; 0 0 1];
        svec_P = term1*term2*term3*[1;0;0];
        se = svec_P[1]; sp = svec_P[2]; sh = svec_P[3]; #breaking it down

        # Ballistic Evolution of state (as a function of xslow)
        e_E = eph.eccentricity
        a_E = eph.semiMajorAxis
        tru_E = get_heliocentric_position(eph, t) # earth heliocentric true anom evaluated at current time
        nudot = (1+e*cos(tru))^2 / (1-e^2)^(3/2) * sqrt(mu/a^3);
        nudot_earth = (1+e_E*cos(tru_E))^2 / (1-e_E^2)^(3/2) * sqrt(mu_sun/a_E^3);

        f0 = @SVector [0; 0; 0; 0; -nudot_earth; nudot];

        # BEST CASE TIME TO GO's
        # Semi-major axis:
        nustar_a = atan(sig_a*se, -sig_a*sp)
        R_H_O_star_a = @SArray [  # Rotation from hill to orbit frame evaluated at nustar
        cos(ape + nustar_a)*cos(lam) - sin(ape + nustar_a)*cos(inc)*sin(lam) cos(ape + nustar_a)*sin(lam) + sin(ape + nustar_a)*cos(inc)*cos(lam) sin(ape + nustar_a)*sin(inc);
         -sin(ape + nustar_a)*cos(lam) - cos(ape + nustar_a)*cos(inc)*sin(lam) cos(ape + nustar_a)*cos(inc)*cos(lam)-sin(ape + nustar_a)*sin(lam) cos(ape + nustar_a)*sin(inc);
         sin(inc)*sin(lam) -cos(lam)*sin(inc) cos(inc)
        ]
        F_a = F(a, e, inc, ape, nustar_a, mu)
        p_a = transpose(-sig_a*eps_a'*F_a*R_H_O_star_a);
        p_ay = p_a[2] 
        p_az = p_a[3]
        betastar_a = atan(-p_ay, -p_az)
        alphastar_a = calculate_alpha_star(p_a, sc)
        alphastar_a = median([params.alpha_min, alphastar_a, params.alpha_max])
        ustar_a = @SVector [alphastar_a; betastar_a]  # elementwise optimal control (EOC)
        astar_hill = aSRP(ustar_a, sc, eph, tru_E) # SRP accel. evaluated at EOC
        adotnn = sig_a*eps_a'*f0 - p_a'*astar_hill
        
        # Ecccentricity
            # For eccentricity, two edotnn's are computed and compared, smaller is taken and used in Q
            # This will be done in the below for loop
        edotnn_set = [0.0, 0.0]
        for n = [0, 1]
            nustar_e = 0.5*atan(sig_e*se, -sig_e*sp) + n*pi
            R_H_O_star_e = hill_to_orbit_transform(inc, ape, lam, nustar_e)
            F_e = F(a, e, inc, ape, nustar_e, mu)
            p_e = transpose(-sig_e*eps_e'*F_e*R_H_O_star_e);
            p_ey = p_e[2]
            p_ez = p_e[3]
            betastar_e = atan(-p_ey, -p_ez)
            alphastar_e = calculate_alpha_star(p_e, sc)
            alphastar_e = median([params.alpha_min, alphastar_e, params.alpha_max])
            ustar_e = @SVector [alphastar_e; betastar_e]  # EOC
            astar_hill = aSRP(ustar_e, sc, eph, tru_E)
            edotnn_set[n+1] = sig_e*eps_e'*f0 - p_e'*astar_hill
        end
        edotnn = min(edotnn_set[1], edotnn_set[2])

        # Inclination
        nustar_inc = pi/2 - ape + sign(sig_inc*sh)*(asin(e*sin(ape))+pi/2)
        return edotnn
end


"""
F: Supplementary function for repeated calls of F(xslow, nustar) as an intermediary step of calculate best-case time-to-go
Inputs: 
a: semimajor axis [km]
e: eccentricity
inc: inclination [rad]
ape: argument of perigee [rad]
tru: trualy [rad]
mu: gravitational parameter 

Outputs:
F: 6x3 Matrix from gauss_variational_eqn
"""
function F(a, e, inc, ape, tru, mu)

    p = a*(1-e^2)
    h = sqrt(mu*p)
    r = p/(1*e*cos(tru))

    out = @SArray [
        2*a^2*e*sin(tru) 2*a^2*p/r 0;
        p*sin(tru) (p+r)*cos(tru)+r*e 0;
        0 0 r*cos(tru+ape);
        -p*cos(tru)/e (p+r)*sin(tru)/e -r*sin(tru+ape)/tan(inc);
        0 0 r*sin(tru+ape)/sin(inc);
        p*cos(tru)/e -(p+r)*sin(tru)/e 0
    ]
    out = out*1/h
    return out
end

"""
calculate_alpha_star: Supplementary function to calculate the term alphastar_0
Solar Sailing Q Law Full Derivation, Equation 20
Inputs: 
p: p_oe vector, calculated at an optimal true anomaly 
"""
function calculate_alpha_star(p, sc::basicSolarSail)
    p_x = p[1] 
    p_y = p[2] 
    p_z = p[3]
    C1 = sc.C[1]
    C2 = sc.C[2]
    C3 = sc.C[3]

    k = p_y/sqrt(p_y^2 + p_z^2);
    alphastar_0 = atan(1/4 * (-3*k + sqrt(8+9*k^2)));

    F_alpha = k*(3*C1*cos(alphastar_0)^2 + 2*C2*cos(alphastar_0) + C3)*sin(alphastar_0) - C1*cos(alphastar_0)*(1-3*sin(alphastar_0)^2) - C2*cos(2*alphastar_0);

    F_alpha_alpha = k*(3*C1*(cos(alphastar_0)-2*cos(alphastar_0))*sin(alphastar_0) + 2*C2*cos(2*alphastar_0) + C3*cos(alphastar_0)) - C1*sin(alphastar_0)*(2-9*cos(alphastar_0))+2*C2*sin(2*alphastar_0);

    alphastar = alphastar_0 - F_alpha/F_alpha_alpha;

    return alphastar
end

"""
aSRP: A functiomn to be called within the Q calculation, made for calculating the SRP acceleration with elementwise optimal controls (and true anomaly)
Notes:
    - SRP acceleration vector is calculated in the sun-centered hill frame
inputs:
    u: control variables [alphastar, betastar]
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
    G0 = 1.02E14 # solar flux constant at earth
    a1 = C1*(cos(α))^2+C2*cos(α)+C3  # sun frame r-direction, not scaled
    a2 = -(C1*cos(α)+C2)*sin(α)sin(β)  # sun frame theta-direction, not scaled
    a3 = -(C1*cos(α)+C2)*sin(α)cos(β)  # sun frame h-direction, not scaled
    a_SRP = cos(α)*sc.areaParam * G0/d^2 * @SVector [a1; a2; a3]  # aSRP expressed in sun frame
    return a_SRP
end

