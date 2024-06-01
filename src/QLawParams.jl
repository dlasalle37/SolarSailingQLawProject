# Type declarations
abstract type QLawType end
abstract type Keplerian <: QLawType end
abstract type Oguri     <: QLawType end 

# Q-Law Parameters Struct
mutable struct QLawParams{T<:QLawType}
    # Basic Parameters
    sc::basicSolarSail
    eph::Ephemeride
    mu::Float64
    mu_sun::Float64

    # Living Parameters
    current_time::Any # time (ephemeris time) Must be able to be a Forward Diff DualType
    alpha::Float64 # control angle alpha [rad]
    beta::Float64 # control angle beta [rad]
    eclipsed::Bool

    #Initial oe's
    oe0::Vector{Float64}

    # Targets 
    oet::Vector{Float64}

    # oe tolerances
    oeTols::Vector{Float64}

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

    # Integration Parameters
    step_size::Float64  # time step size 
    max_sim_time::Float64 # max # of seconds to simulate for
    abstol::Float64
    reltol::Float64

    # Writing data
    writeData::Bool

end


# Constructor 
function QLawParams(
    sc::basicSolarSail,
    eph::Ephemeride,
    oe0,
    oet,
    oeTols; ### MANDATORY INPUTS END HERE ###
    mu::Float64=GRAV_PARAMS[399], #Defaults to earth
    mu_sun::Float64=GRAV_PARAMS[10], #helpful to carry this around, not totally necessary
    Woe::Vector{Float64}=[1.0; 1.0; 1.0; 1.0; 1.0],
    Wp::Int=1,
    Aimp::Int=1,
    kimp::Int=100,
    Aesc::Int=1,
    kesc::Int=10,
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
    alpha=0.0,
    step_size=60.0, #seconds
    max_sim_time = 100*86400.0, # seconds (days*sec/day)
    abstol=1.0E-6,
    reltol=1.0E-6,
    writeData=true,
    type=Oguri
)
    current_time = eph.t0

    # Find out if initial conditions are eclipsed
    coee = getCOE(eph, current_time)
    nue = coee[6]
    (r, ~) = coe2rv(oe0[1], oe0[2], oe0[3], oe0[4], oe0[5]+nue, oe0[6], mu)
    eclipsed = isEclipsed(eph, r, current_time)

    qlawparam = QLawParams{type}(sc, eph, mu, mu_sun, current_time, alpha, beta, eclipsed, oe0, oet, oeTols, Woe, Wp, 
        Aimp, kimp, Aesc, kesc, rp_min, a_esc, m_petro, n_petro, r_petro, alpha_min, alpha_max, beta_min, beta_max, step_size, max_sim_time, abstol, reltol,
        writeData)

    return qlawparam
end