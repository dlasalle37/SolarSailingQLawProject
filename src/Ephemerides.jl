"""
    Ephemeride

Ephemeride struct for storing ephemeris data for a single target body.

# Fields
- `t0::Float64`: Initial epoch of ephemeride in SPICE ET (TBD)
- `tf::Float64`: Final epoch of ephemeride in SPICE ET (TBD)
- `targ::Int`: Target NAIF ID
- `mu::Float64`: Target gravitational parameter
- `obs::Int`: Observer NAIF ID
- `frame::String`: Reference frame string (e.g. "J2000")
- `abcorr::String`: Aberation corrections (e.g. "NONE")
- `spline::CubicSpline`: Cubic spline interpolant
"""
struct Ephemeride
    # Initial and final epochs of ephemeride in SPICE ET (TBD)
    t0::Float64
    tf::Float64

    # Target NAIF ID
    targ::Int

    # Target gravitational parameter
    mu::Float64

    # Observer NAIF ID
    obs::Int

    # Reference frame string
    frame::String

    # Aberation corrections
    abcorr::String

    # Cubic spline interpolant
    spline::CubicSpline
end

"""
    Ephemerides

Ephemerides struct for storing ephemeris data for multiple target bodies.

# Fields
- `t0::Float64`: Initial epoch of ephemerides in SPICE ET (TBD)
- `tf::Float64`: Final epoch of ephemerides in SPICE ET (TBD)
- `targIDs::Vector{Int}`: Vector of target NAIF IDs
- `obs::Int`: Observer NAIF ID
- `ephems::Vector{Ephemeride}`: Vector of Ephemeride structs
- `cache::Vector{Vector{Float64}}`: Cache for interpolated states
"""
struct Ephemerides
    # Initial and final epochs in ephemeride in SPICE ET (TBD)
    t0::Float64
    tf::Float64

    # Target IDs
    targIDs::Vector{Int}

    # Observer ID
    obs::Int

    # Vector of Ephemerid structs
    ephems::Vector{Ephemeride}

    # Cache for interpolated states
    cache::Vector{Vector{Float64}}
end

"""
    Ephemerides(tspan::Tuple, n_points::Int, targs, obs, ref; abcorr = "NONE")

Constructor for Ephemerides struct. Computes ephemerides for multiple target

# Arguments
- `tspan::Tuple`: Tuple of initial and final epochs in SPICE ET (TBD)
- `n_points::Int`: Number of points to employ in interpolants
- `targs`: Vector of target NAIF IDs
- `obs`: Observer NAIF ID
- `ref`: Reference frame string (e.g. "J2000")

# Keywords
- `abcorr::String="NONE"`: Aberation corrections

# Returns
- `ephems::Ephemerides`: Ephemerides struct
"""
function Ephemerides(tspan::Tuple, n_points::Int, targs, obs, ref; abcorr = "NONE")
    # Create vector of Ephemerides
    ephems = [Ephemeride(
        tspan, n_points, targs[i], obs, ref; abcorr = abcorr,
    ) for i in eachindex(targs)]

    # Create cache
    cache = [zeros(3) for i in eachindex(targs)]

    # Return constructed Ephemerides
    return Ephemerides(tspan[1], tspan[2], targs, obs, ephems, cache)
end

"""
    Ephemerides(tspan::Tuple, n_points::Int, targs, obs, ref; abcorr = "NONE")

Constructor for Ephemerides struct. Computes ephemerides for multiple target using a
    different number of points for each target.

# Arguments
- `tspan::Tuple`: Tuple of initial and final epochs in SPICE ET (TBD)
- `n_points::AbstractVector{Int}`: Vecter of integers indicating the number of points to
    employ for each interpolant
- `targs`: Vector of target NAIF IDs
- `obs`: Observer NAIF ID
- `ref`: Reference frame string (e.g. "J2000")

# Keywords
- `abcorr::String="NONE"`: Aberation corrections

# Returns
- `ephems::Ephemerides`: Ephemerides struct

# Throws
- `ArgumentError`: If length of `n_points` is not equal to length of `targs`
"""
function Ephemerides(tspan::Tuple, n_points::AbstractVector{Int}, targs, obs, ref; abcorr = "NONE")
    # Check that length n_points is equal to length targs
    if length(n_points) != length(targs)
        throw(ArgumentError("If setting n_points for each target body, " *
            "n_points must be the same length as targs."))
    end

    # Create vector of Ephemerides
    ephems = [Ephemeride(
        tspan, n_points[i], targs[i], obs, ref; abcorr = abcorr,
    ) for i in eachindex(targs)]

    # Create cache
    cache = [zeros(3) for i in eachindex(targs)]

    # Return constructed Ephemerides
    return Ephemerides(tspan[1], tspan[2], targs, obs, ephems, cache)
end

"""
    Ephemeride(tspan::Tuple, n_points, targ, obs, ref; abcorr = "NONE")

Constructor for Ephemeride struct. Computes ephemeride for a single target body.

# Arguments
- `tspan::Tuple`: Tuple of initial and final epochs in SPICE ET (TBD)
- `n_points::Int`: Number of points to employ in interpolant
- `targ`: Target NAIF ID
- `obs`: Observer NAIF ID
- `ref`: Reference frame string (e.g. "J2000")

# Keywords
- `abcorr::String="NONE"`: Aberation corrections

# Returns
- `ephem::Ephemeride`: Ephemeride struct
"""
function Ephemeride(tspan::Tuple, n_points, targ, obs, ref; abcorr = "NONE")
    # Compute time stamps for interpolants
    ts = LinRange(0.0, 1.0, n_points)

    # Get body gravitational parameter
    mu = get_gm(targ)

    # Get bodies position at each epoch
    states = zeros(n_points, 3)
    for i in 1:n_points
        pos, ld = spkezp(targ, tspan[1] + (tspan[2] - tspan[1])*ts[i], ref, abcorr, obs)
        states[i,:] .= pos
    end

    # Get bodies velocity at initial and final epoch
    state0, ld = spkez(targ, tspan[1], ref, abcorr, obs)
    statef, ld = spkez(targ, tspan[2], ref, abcorr, obs)
    dx0dt      = view(state0, 4:6)
    dxfdt      = view(statef, 4:6)

    # Convert to partials wrt scaled time
    sf         = tspan[2] - tspan[1]
    dx0dτ      = SVector(sf*dx0dt[1], sf*dx0dt[2], sf*dx0dt[3])
    dxfdτ      = SVector(sf*dxfdt[1], sf*dxfdt[2], sf*dxfdt[3])

    # Compute cubic spline interpolant
    spline = CubicSpline(ts, states, dx0dτ, dxfdτ)

    return Ephemeride(tspan[1], tspan[2], targ, mu, obs, ref, abcorr, spline)
end

"""
    get_gravitational_parameter(e::Ephemeride)

Returns the gravitational parameter of the target body.

# Arguments
- `e::Ephemeride`: Ephemeride struct

# Returns
- `mu::Float64`: Gravitational parameter of target body
"""
function get_gravitational_parameter(e::Ephemeride)
    return e.mu
end

"""
    get_gravitational_parameter(e::Ephemerides)

Returns the gravitational parameter of each target body.

# Arguments
- `e::Ephemerides`: Ephemerides struct

# Returns
- `mus::NTuple{N,Float64}`: Tuple of gravitational parameters for each target body
"""
function get_gravitational_parameter(e::Ephemerides)
    return ntuple(i -> get_gravitational_parameter(e.ephems[i]), length(e.ephems))
end

"""
    interpolate(e::Ephemeride, t)

Interpolates the ephemeride to time t.

# Arguments
- `e::Ephemeride`: Ephemeride struct
- `t::Float64`: Time to interpolate to in SPICE ET (TBD)

# Returns
- `state::Vector{Float64}`: Interpolated state

# Throws
- `ArgumentError`: If t is outside of the Ephemeride's domain
"""
function interpolate(e::Ephemeride, t)
    # Check that t is in bounds
    if t < e.t0 || t > e.tf
        throw(ArgumentError("t is outside of the Ephemeride's domain."))
    end

    # Compute scaled time
    τ = (t - e.t0) / (e.tf - e.t0)

    # Return interpolant
    return interpolate(e.spline, τ)
end

"""
    interpolate(e::Ephemerides, t)

Interpolates the ephemerides to time t.

# Arguments
- `e::Ephemerides`: Ephemerides struct
- `t::Float64`: Time to interpolate to in SPICE ET (TBD)

# Returns
- `states::Vector{Vector{Float64}}`: Interpolated states
"""
function interpolate(e::Ephemerides, t)
    # Check that t is in bounds
    if t < e.t0 || t > e.tf
        throw(ArgumentError("t is outside of the Ephemeride's domain."))
    end

    # Compute scaled time
    τ = (t - e.t0) / (e.tf - e.t0)

    # Grab cache for easy usage
    cache = e.cache

    # Interpolate for each body
    @inbounds for i in eachindex(cache)
        cache[i] .= interpolate(e.ephems[i], t)
    end
    return cache
end

"""
    interpolate(e::Ephemerides, targ, t)

Interpolates the ephemerides to time t.

# Arguments
- `e::Ephemerides`: Ephemerides struct
- `targ::Int`: Target NAIF ID
- `t::Float64`: Time to interpolate to in SPICE ET (TBD)

# Returns
- `states::SVector`: Interpolated states
"""
function interpolate(e::Ephemerides, targ, t)
    # Check that targ is included in ephemerides
    if !(targ in e.targIDs)
        throw(ArgumentError("Ephemerides does not contain ephemeris for target ID: $targ"))
    end

    # Return interpolated state
    for i in eachindex(e.ephems)
        if e.ephems[i].targ == targ
            return interpolate(e.ephems[i], t)
        end
    end
end

"""
    create_serialized_ephemerides(utc0, utcf, n_points, targs, obs, ref)

Creates and saves ephemerides in a serialized form. This is used to avoid SPICE calls
    during maneuver plan computation. The serialized ephemerides are saved at
    `.../DSMManeuverPlanner/data/ephem/ephem.jld2`.

# Arguments
- `utc0::String`: Initial epoch in UTC
- `utcf::String`: Final epoch in UTC
- `n_points::Int`: Number of points to employ in interpolants
- `targs`: Vector of target NAIF IDs
- `obs`: Observer NAIF ID
- `ref`: Reference frame string (e.g. "J2000")

# Returns
- `nothing`

# Examples
```jldoctest
julia> utc0    = "2028 January 1, 00:00:00.0";
julia> utcf    = "2029 March 1, 00:00:00.0";
julia> obs     = 602;
julia> targs   = [699];
julia> npoints = 10000;
julia> create_serialized_ephemerides(utc0, utcf, npoints, targs, obs, "J2000");
```
"""
function create_serialized_ephemerides(utc0, utcf, n_points, targs, obs, ref)
    # Load spice kernels
    furnsh_kernals()

    # Compute time stams in ET
    et0 = str2et(utc0)
    etf = str2et(utcf)

    # Create Ephemerides
    ephem = Ephemerides((et0,etf), n_points, targs, obs, ref)

    # Unload kernals
    unload_kernals()

    # Save with JLD2
    path = joinpath(@__DIR__, "..", "data", "ephem", "ephem.jld2") 
    save_object(path, ephem)
    return ephem
end

"""
    load_serialized_ephemerides()

Loads serialized ephemerides from `../data/ephem/ephem.jld2`.

# Returns
- `ephem::Ephemerides`: Ephemerides struct

# Examples
```jldoctest
julia> ephem = load_serialized_ephemerides();
```
"""
function load_serialized_ephemerides()
    path = joinpath(@__DIR__, "..", "data", "ephem", "ephem.jld2")
    return load_object(path)
end

"""
    getState
Use the interpolation tools to retrieve the entire state of a target body at time some time t.

# Inputs:
    - e::Ephemeride; instantiated ephemeride struct
    - t: some time in SPICE Ephemeris time between e.t0 and e.Transfer
# Outputs:
    - state: 6x1 state vector [pos{3x1}; vel{3x1}] at time t with units [km{3x1};km/s{3x1}]
"""
function getState(e::Ephemeride, t)
    # Check bounds of time
    if t < eph.t0 || t > e.tf
        throw(DomainError("Time value given to getState() is out of the Ephemeride's bounds"))
    end

    spline = e.spline
    # Scale time
    τ = (t - e.t0) / (e.tf - e.t0) 

    r = interpolate(spline, τ)

    drdτ = getPositionPartials(spline, τ)
    dτdt = 1/(e.tf-e.t0)
    v = drdτ*dτdt 

    return SVector([r;v])

end