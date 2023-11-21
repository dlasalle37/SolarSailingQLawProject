"""
compute_control: function to compute alphastar and betastar at a given instant in time
inputs: 
    -x: state vector [a, e, inc, ape, lam, tru] w/ units [km, none, rad, rad, rad, rad]
    -dQdx: gradient of Q wrt state vector
outputs: 
    -alphastar: control variable alpha at given time instant
    -betastar: control variable beta at given time instant
"""
function compute_control(x, params::QLawParams)
    # Unpacking:
    mu = params.mu
    a = x[1]
    e = x[2]
    inc = x[3]
    ape = x[4]
    lam = x[5]
    tru = x[6]

    sc=params.sc
    C1 = sc.C[1]
    C2 = sc.C[2]
    C3 = sc.C[3]

    
    # Keep inc, e away from zero
    if e<= 1.0E-4
        e = 1.0E-4
    end
    if inc <= 1.0E-4
        inc = 1.0E-4
    end

    # Calculation
    dQdx = FiniteDiff.finite_difference_gradient(x->calculate_Q(x, params), x)
    Fx = F(a, e, inc, ape, tru, mu)
    R_H_O = hill_to_orbit_transform(inc, ape, lam, tru)  # rotation matrix for current states
    pvec = -transpose(dQdx'*Fx*R_H_O)
    px = pvec[1]
    py = pvec[2]
    pz = pvec[3]
    
    # Alphastar
    k = px/sqrt(py^2+pz^2)
    αstar0 = atan(0.25*(-3*k+sqrt(8+9*k^2)))
    Fstar = k*(3*C1*cos(αstar0)^2+2*C2*cos(αstar0)+C3)*sin(αstar0)-C1*cos(αstar0)*(1-3*sin(αstar0)^2)-C2*cos(2*αstar0)
    Fstar_prime = k*(3*C1*(cos(αstar0)-2*cos(αstar0)*sin(αstar0)) + 2*C2*cos(2*αstar0)+C3*cos(αstar0)) - 
        C1*sin(αstar0)*(2-9*cos(αstar0)) + 2*C2*sin(2*αstar0)
    alphastar = αstar0 - Fstar/Fstar_prime
    alphastar = median([params.alpha_min, alphastar, params.alpha_max]) # enforcing alphastar range constraint

    # betastar
    betastar = atan(-py, -pz)
    return alphastar, betastar, dQdx
end

"""
The below function is where the actual calculation of Q is done
"""
function calculate_Q(x, params)

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
        oet = params.oet
        a_t = oet[1]
        e_t = oet[2]
        inc_t = oet[3]
        ape_t = oet[4]
        lam_t = oet[5]

        # Weights
        Woe = params.Woe
        Wa = Woe[1]
        We = Woe[2]
        Winc = Woe[3]
        Wape = Woe[4]
        Wlam = Woe[5]

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
        P = P1 + P2

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
        svec_P = @SArray [
            cos(lam)*cos(ape)-sin(lam)*cos(inc)*sin(ape);
            -cos(lam)*sin(ape)-sin(lam)*cos(inc)*cos(ape);
            sin(lam)*sin(inc)
        ]
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
        R_H_O_star_a = hill_to_orbit_transform(inc, ape, lam, nustar_a)
        F_a = F(a, e, inc, ape, nustar_a, mu)
        p_a = -transpose(sig_a*eps_a'*F_a*R_H_O_star_a);
        p_ay = p_a[2] 
        p_az = p_a[3]
        betastar_a = atan(-p_ay, -p_az)
        alphastar_a = calculate_alpha_star(p_a, sc)
        alphastar_a = median([params.alpha_min, alphastar_a, params.alpha_max]) # enforcing alphastar range constraint
        ustar_a = @SVector [alphastar_a; betastar_a]  # elementwise optimal control (EOC)
        astar_hill = aSRP(ustar_a, sc, eph, tru_E) # SRP accel. evaluated at EOC
        adotnn = sig_a*eps_a'*f0 - p_a'*astar_hill # positive denominator of best-case ttg for A
        tau_a = abs(dista)/-adotnn  # best-case ttg term
        
        # Ecccentricity
            # For eccentricity, two edotnn's are computed and compared, smaller is taken and used in Q
            # This will be done in the below for loop
        edotnn_set = [0.0, 0.0]
        for n = [0, 1]
            nustar_e = 0.5*atan(sig_e*se, -sig_e*sp) + n*pi
            R_H_O_star_e = hill_to_orbit_transform(inc, ape, lam, nustar_e)
            F_e = F(a, e, inc, ape, nustar_e, mu)
            p_e = -transpose(sig_e*eps_e'*F_e*R_H_O_star_e);
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
        tau_e = abs(diste)/-edotnn

        # Inclination
        nustar_inc = pi/2 - ape + sign(sig_inc*sh)*(asin(e*sin(ape))+pi/2)
        R_H_O_star_inc = hill_to_orbit_transform(inc, ape, lam, nustar_inc)
        F_inc = F(a, e, inc, ape, nustar_inc, mu)
        p_inc = -transpose(sig_inc*eps_inc'*F_inc*R_H_O_star_inc);
        p_incy = p_inc[2]
        p_incz = p_inc[3]
        betastar_inc = atan(-p_incy, -p_incz)
        alphastar_inc = calculate_alpha_star(p_inc, sc)
        alphastar_inc = median([params.alpha_min, alphastar_inc, params.alpha_max])
        ustar_inc = @SVector [alphastar_inc; betastar_inc]
        astar_hill = aSRP(ustar_inc, sc, eph, tru_E)
        incdotnn = sig_inc*eps_inc'*f0 - p_inc'*astar_hill
        tau_inc = abs(distinc)/-incdotnn

        # Argument of periapsis
        if Wape == 0.0
            tau_ape = 0.0
        else # only calculate this if Wape != 0
            # This one is done with a brute force method:
            nu = range(0, 2*pi, 200)  # 200 evenly spaced points from 0 to 2*pi
            mappedvals = @MArray zeros(length(nu)) # initialize solution
            for i in range(1, length(nu))
                mappedvals[i] = sig_ape*p/(h*(1+e*cos(nu[i]))) * (1/e * sin(nu[i])*(-se*sin(nu[i])+sp*cos(nu[i])) - sin(nu[i]+ape)/tan(inc) * sh)
            end
            nustar_ape = nu[argmin(mappedvals)] # the nu that minimized mappedvals is the approx. nustar for ape

            # From here, it is the same as the others
            R_H_O_star_ape = hill_to_orbit_transform(inc, ape, lam, nustar_ape)
            F_ape = F(a, e, inc, ape, nustar_ape, mu)
            p_ape = -transpose(sig_ape*eps_ape'*F_ape*R_H_O_star_ape);
            p_apey = p_ape[2]
            p_apez = p_ape[3]
            betastar_ape = atan(-p_apey, -p_apez)
            alphastar_ape = calculate_alpha_star(p_ape, sc)
            alphastar_ape = median([params.alpha_min, alphastar_ape, params.alpha_max])
            ustar_ape = @SVector [alphastar_ape; betastar_ape]
            astar_hill = aSRP(ustar_ape, sc, eph, tru_E)
            apedotnn = sig_ape*eps_ape'*f0 - p_ape'*astar_hill
            tau_ape = abs(distape)/-apedotnn
        end

        # Longitude of ascending node from ir (lambda)
        if Wlam == 0.0
            tau_lam = 0.0
        else # only calculate if Wlam !=0
            nustar_lam = pi - ape + sign(sig_lam*sh/sin(inc))*acos(e*cos(ape))
            R_H_O_star_lam = hill_to_orbit_transform(inc, ape, lam, nustar_lam)
            F_lam = F(a, e, inc, ape, nustar_lam, mu)
            p_lam = -transpose(sig_lam*eps_lam'*F_lam*R_H_O_star_lam);
            p_lamy = p_lam[2]
            p_lamz = p_lam[3]
            betastar_lam = atan(-p_lamy, -p_lamz)
            alphastar_lam = calculate_alpha_star(p_lam, sc)
            alphastar_lam = median([params.alpha_min, alphastar_lam, params.alpha_max])
            ustar_lam = @SVector [alphastar_lam; betastar_lam]
            astar_hill = aSRP(ustar_lam, sc, eph, tru_E)
            lamdotnn = sig_lam*eps_lam'*f0 - p_lam'*astar_hill
            tau_lam = abs(distlam)/-lamdotnn
        end

        ###################################################################################################################################################
        # PUTTING Q TOGETHER :)
        Q = (1+Wp*P)*(Wa*Sa*tau_a^2 + We*tau_e^2 + Winc*tau_inc^2 + Wape*tau_ape^2 + Wlam*tau_lam^2)
        ###################################################################################################################################################

        return Q
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

    k = p_x/sqrt(p_y^2 + p_z^2);
    alphastar_0 = atan(1/4 * (-3*k + sqrt(8+9*k^2)));

    F_alpha = k*(3*C1*cos(alphastar_0)^2 + 2*C2*cos(alphastar_0) + C3)*sin(alphastar_0) - C1*cos(alphastar_0)*(1-3*sin(alphastar_0)^2) - C2*cos(2*alphastar_0);

    F_alpha_alpha = k*(3*C1*(cos(alphastar_0)-2*cos(alphastar_0)*sin(alphastar_0)) + 2*C2*cos(2*alphastar_0) + C3*cos(alphastar_0)) - C1*sin(alphastar_0)*(2-9*cos(alphastar_0))+2*C2*sin(2*alphastar_0);

    alphastar = alphastar_0 - F_alpha/F_alpha_alpha;

    return alphastar
end


"""
This is a function for debugging
    analytical partial for penalty wrt eccentricity
"""
function analytical_dPde(e, a, params::QLawParams)
    kimp = params.kimp
    Aimp = params.Aimp
    rpmin = params.rp_min
    dPde = kimp*a/rpmin * Aimp*exp(kimp * (1-a*(1-e)/rpmin))

end

"""
use finitediff on this and compare to above
"""
function Pfun(e, a, params)
    kimp = params.kimp
    Aimp = params.Aimp
    rpmin = params.rp_min
    P = Aimp*exp(kimp * (1-a*(1-e)/rpmin))
end