mutable struct QLawParams
    # Basic Parameters
    sc::basicSolarSail
    eph::TwoBodyEphemeride
    mu::Float64
    mu_sun::Float64

    # Living Parameters
    current_time::Float64 # time (seconds past j2000) (this should be continually updated through zero-order hold)
    alpha::Float64 # control angle alpha [rad]
    beta::Float64 # control angle beta [rad]

    #Initial oe's
    oe0::Vector{Float64}

    # Targets 
    oet::Vector{Float64}

    # oe tolerances
    oetol::Vector{Float64}

    # Weights
    Woe::Vector{Float64}

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
    oe0,
    oet,
    oetol; ### MANDATORY INPUTS END HERE ###
    mu::Float64=EARTH_MU,
    mu_sun::Float64=SUN_MU,
    Woe::Vector{Float64}=[1.0; 1.0; 1.0; 1.0; 1.0],
    Wp::Int=1,
    Aimp::Int=1,
    kimp::Int=100,
    Aesc::Int=1,
    kesc::Int=100,
    rp_min::Float64=6578.0,
    a_esc::Float64=1.0E5,
    m_petro::Int=3,
    n_petro::Int=4,
    r_petro::Int=2,
    alpha_min::Float64=0.0,
    alpha_max::Float64=pi/2,
    beta_min::Float64=-pi,
    beta_max::Float64=1.0*pi,
    beta=0.0,
    alpha=0.0
)

    qlawparam = QLawParams(sc, eph, mu, mu_sun, current_time, alpha, beta, oe0, oet, oetol, Woe, Wp, 
        Aimp, kimp, Aesc, kesc, rp_min, a_esc, m_petro, n_petro, r_petro, alpha_min, alpha_max, beta_min, beta_max)

    return qlawparam
end