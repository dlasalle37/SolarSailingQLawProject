mutable struct QLawParams
    # Basic Parameters
    sc::basicSolarSail
    eph::TwoBodyEphemeride
    mu::Float64
    mu_sun::Float64
    current_time::Float64 # time (seconds past j2000)
    # Targets 
    a_t::Float64
    e_t::Float64
    inc_t::Float64
    ape_t::Float64
    lam_t::Float64

    # Weights
    Wa::Int
    We::Int
    Winc::Int
    Wape::Int
    Wlam::Int

    # Penalty Parameters
    Wp::Int # Penalty weighting terms
    Aimp::Int # magnitude, impact_constraint
    kimp::Int # sharpness, impact_constraint
    Aesc::Int # magnitude, escape_constraint
    kesc::Int # sharpness, escape_constraint
    rp_min::Float64
    a_esc::Float64

    # Scaling
    m_petro::Int
    n_petro::Int
    r_petro::Int

    # Control Bounds
    alpha_min::Float64
    alpha_max::Float64
    beta_min::Float64
    beta_max::Float64
end


# Constructor 
function QLawParams(
    sc::basicSolarSail,
    eph::TwoBodyEphemeride,
    current_time::Float64,
    a_t,
    e_t,
    inc_t,
    ape_t,
    lam_t; ### MANDATORY INPUTS END HERE ###
    mu::Float64=EARTH_MU,
    mu_sun::Float64=SUN_MU,
    Wa::Int=1,
    We::Int=1,
    Winc::Int=1,
    Wape::Int=1,
    Wlam::Int=1,
    Wp::Int=1,
    Aimp::Int=1,
    kimp::Int=100,
    Aesc::Int=1,
    kesc::Int=100,
    rp_min::Float64=6378.0,
    a_esc::Float64=1.0E5,
    m_petro::Int=3,
    n_petro::Int=4,
    r_petro::Int=2,
    alpha_min::Float64=0.0,
    alpha_max::Float64=pi/2,
    beta_min::Float64=-pi,
    beta_max::Float64=1.0*pi
)
    qlawparam = QLawParams(sc, eph, mu, mu_sun, current_time, a_t, e_t, inc_t, ape_t, lam_t, Wa, We, Winc, Wape, Wlam, Wp, 
        Aimp, kimp, Aesc, kesc, rp_min, a_esc, m_petro, n_petro, r_petro, alpha_min, alpha_max, beta_min, beta_max)

    return qlawparam
end