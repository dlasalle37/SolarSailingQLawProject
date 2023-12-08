"""
compute_control: function to compute alphastar and betastar at a given instant in time
inputs: 
    -x: state vector [a, e, inc, ape, lam, tru] w/ units [km, none, rad, rad, rad, rad]
    -params: QLaw Params struct containing all supplementary info
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

    # DQDX Calculation:
    
    #= 
    There are two different numerical differentiation methods in the 3-line block below
      To use ForwardDiff: uncomment first two lines, comment third line
      To use FiniteDiff: uncomment third line, comment first two lines  
    =#
    cfg = ForwardDiff.GradientConfig(calculate_Q, x) # GradientConfig
    dQdx = ForwardDiff.gradient(x->calculate_Q(x, params), x, cfg, Val{false}()) 
    #dQdx = FiniteDiff.finite_difference_gradient(x->calculate_Q(x, params), x) # using finite diff
    
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
    Fstar_prime = k*(3*C1*(cos(αstar0)^3-2*cos(αstar0)*sin(αstar0)^2) + 2*C2*cos(2*αstar0)+C3*cos(αstar0)) - 
        C1*sin(αstar0)*(2-9*cos(αstar0)^2) + 2*C2*sin(2*αstar0)
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
        if mu*p <0
            @infiltrate
        end
        h = sqrt(mu*p)
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
        svec_P = @SArray [
            cos(lam)*cos(ape)-sin(lam)*cos(inc)*sin(ape);
            -cos(lam)*sin(ape)-sin(lam)*cos(inc)*cos(ape);
            sin(lam)*sin(inc)
        ]
        se = svec_P[1]; sp = svec_P[2]; sh = svec_P[3]; #breaking it down

        # Ballistic Evolution of state (as a function of x)
        e_E = eph.eccentricity
        a_E = eph.semiMajorAxis
        tru_E = get_heliocentric_position(eph, t) # earth heliocentric true anom evaluated at current time
        nudot = (1+e*cos(tru))^2 / (1-e^2)^(3/2) * sqrt(mu/a^3);
        nudot_earth = (1+e_E*cos(tru_E))^2 / (1-e_E^2)^(3/2) * sqrt(mu_sun/a_E^3);

        f0 = @SVector [0; 0; 0; 0; -nudot_earth; nudot];

        # BEST CASE TIME TO GO's
        # Semi-major axis:
        nustar_a = atan(sig_a*se, -sig_a*sp)
        adotnn = oedotnn(a, e, inc, ape, lam, tru, nustar_a, sig_a, eps_a, tru_E, f0, params)
        tau_a = abs(dista)/-adotnn  # best-case ttg term
        
        # Ecccentricity
            # For eccentricity, two edotnn's are computed and compared, smaller is taken and used in Q
        nustar_e1 = 0.5*atan(sig_e*se, -sig_e*sp) + 0*pi
        nustar_e2 = 0.5*atan(sig_e*se, -sig_e*sp) + 1*pi
        edotnn1 = oedotnn(a, e, inc, ape, lam, tru, nustar_e1, sig_e, eps_e, tru_E, f0, params)
        edotnn2 = oedotnn(a, e, inc, ape, lam, tru, nustar_e2, sig_e, eps_e, tru_E, f0, params)
        edotnn = min(edotnn1, edotnn2)

        # need to clip edotnn if it is near-zero 
        d = distance_to_sun(eph, tru_E)
        G0 = 1.02E14 # solar flux constant at earth [kgkm/s^2]
        a_over_m = sc.areaParam
        ε = 1.0E-5 # given in the Oguri paper as a constant for numerical simulations
        clip_val = -ε*a_over_m*(G0/d^2)*sqrt(a_t/mu)
        if abs(edotnn) <= abs(clip_val) # edotnn should be negative, and needs to be clipped if nearer to zero than clip_val
            edotnn = clip_val
        end

        tau_e = abs(diste)/-edotnn
        
        # Inclination
        nustar_inc = pi/2 - ape + sign(sig_inc*sh)*(asin(e*sin(ape))+pi/2)
        incdotnn = oedotnn(a, e, inc, ape, lam, tru, nustar_inc, sig_inc, eps_inc, tru_E, f0, params)
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
            apedotnn = oedotnn(a, e, inc, ape, lam, tru, nustar_ape, sig_ape, eps_ape, tru_E, f0, params)
            tau_ape = abs(distape)/-apedotnn
        end

        # Longitude of ascending node from ir (lambda)
        if Wlam == 0.0
            tau_lam = 0.0
        else # only calculate if Wlam !=0
            nustar_lam = pi - ape + sign(sig_lam*sh/sin(inc))*acos(e*cos(ape))
            lamdotnn = oedotnn(a, e, inc, ape, lam, tru, nustar_lam, sig_lam, eps_lam, tru_E, f0, params)
            tau_lam = abs(distlam)/-lamdotnn
        end

        ###################################################################################################################################################
        # PUTTING Q TOGETHER :)
        Q = (1+Wp*P)*(Wa*Sa*tau_a^2 + We*tau_e^2 + Winc*tau_inc^2 + Wape*tau_ape^2 + Wlam*tau_lam^2)
        ###################################################################################################################################################
        @infiltrate true
        return Q
end

"""
oedotnn: calculate the oedotnn term in best-case time-to-go based on nustar_oe
INPUTS:
    a: semi-major axis[km]
    e: eccentricity
    inc: inclination [rad]
    ape: arg. periapsis [rad]
    lam: longitude of ascending node [rad]
    nustar_oe: optimal true anomaly for element in oe
    sig_oe: sigma function for element in oe
    eps_oe: selection vector for given element in oe (e.g. [0 1 0 0 0 0] for e)
    tru_E: earth true anomaly at time of calculation [rad]
    f0: ballistic evolution vector at time of calculation [rad/s]
    params: QLawParams struct containing supplementary info
OUTPUT:
    oedotnn: positive denominator of best-case time-to-go for given element in oe
"""
function oedotnn(a, e, inc, ape, lam, tru, nustar_oe, sig_oe, eps_oe, tru_E, f0, params::QLawParams)
    sc = params.sc
    mu = params.mu
    mu_sun = params.mu_sun
    eph = params.eph
    R_H_O_star = hill_to_orbit_transform(inc, ape, lam, nustar_oe)
    Foe = F(a, e, inc, ape, nustar_oe, mu)
    pvec = -transpose(sig_oe*eps_oe'*Foe*R_H_O_star);
    py = pvec[2] 
    pz = pvec[3]
    betastar = atan(-py, -pz)
    alphastar = calculate_alpha_star(pvec, sc)
    alphastar = median([params.alpha_min, alphastar, params.alpha_max]) # enforcing alphastar range constraint
    ustar = @SVector [alphastar; betastar]  # elementwise optimal control (EOC)
    astar_hill = aSRP(ustar, sc, eph, tru_E) # SRP accel. evaluated at EOC
    oedotnn = sig_oe*eps_oe'*f0 - pvec'*astar_hill # positive denominator of best-case ttg for A
    return oedotnn
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
function F(a, e, i, ω, θ, mu)
    p = a*(1-e^2)
    h = sqrt(mu*p)
    r = p/(1+e*cos(θ))

    F = 
    @SArray [
        2*a^2*e*sin(θ) 2*a^2*p/r 0;
        p*sin(θ) (p+r)*cos(θ)+r*e 0;
        0 0 r*cos(θ+ω);
        -p*cos(θ)/e (p+r)*sin(θ)/e -r*sin(θ+ω)/tan(i);
        0 0 r*sin(θ+ω)/sin(i);
        p*cos(θ)/e -(p+r)*sin(θ)/e 0;
    ]
    out = F*1/h
    return out
end

"""
calculate_alpha_star: Supplementary function to calculate the term alphastar_0
Solar Sailing Q Law Full Derivation, Equation 20
Inputs: 
p: p_oe vector, calculated at an optimal true anomaly 
"""
function calculate_alpha_star(p, sc::basicSolarSail)
    px = p[1] 
    py = p[2] 
    pz = p[3]
    C1 = sc.C[1]
    C2 = sc.C[2]
    C3 = sc.C[3]

    k = px/sqrt(py^2+pz^2)
    αstar0 = atan(0.25*(-3*k+sqrt(8+9*k^2)))
    Fstar = k*(3*C1*cos(αstar0)^2+2*C2*cos(αstar0)+C3)*sin(αstar0)-C1*cos(αstar0)*(1-3*sin(αstar0)^2)-C2*cos(2*αstar0)
    Fstar_prime = k*(3*C1*(cos(αstar0)^3-2*cos(αstar0)*sin(αstar0)^2) + 2*C2*cos(2*αstar0)+C3*cos(αstar0)) - 
        C1*sin(αstar0)*(2-9*cos(αstar0)^2) + 2*C2*sin(2*αstar0)
    alphastar = αstar0 - Fstar/Fstar_prime

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