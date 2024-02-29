#=
COMPUTING ERROR BETWEEN SPICE CALLS AND USING INTERPOLATION FOR PLANET POSITIONS
=#

using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

import GLMakie as GM

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP


date = "2023-03-01T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-05T12:30:00"
endEpoch = utc2et(endDate)  # end date


# Create the interpolated Ephemeride
tspan = (epoch, endEpoch)
nInterpolants = 1000
targ = 399
obs = 10
frame = "J2000"
eph = Ephemeride(tspan, nInterpolants, targ, obs, frame)

# Create a vector of interpolated points
nPoints = 10*nInterpolants
r_interp = zeros(nPoints, 3)
t = LinRange(tspan[1], tspan[2], nPoints)
for row in axes(r_interp, 1)
    r_interp[row,:] .= interpolate(eph, t[row])
end

# Call Spice to retrieve position at each point
r_spice = zeros(nPoints, 3)
for row in axes(r_spice, 1)
    tnow = t[row]
    pos, lt = spkezp(targ, tnow, frame, "none", obs)
    r_spice[row,:] .= pos
end

# Compute error
errX = r_interp[:,1] - r_spice[:,1]
errY = r_interp[:,2] - r_spice[:,2]
errZ = r_interp[:,3] - r_spice[:,3]

# Plot error
fig = GM.Figure()
axx = GM.Axis(
    fig[1,1];
    xlabel = "Time[s]",
    ylabel = "Error X [km]"
)
x = GM.lines!(axx, t.-tspan[1], errX)
axy = GM.Axis(
    fig[2,1],
    xlabel="Time[s]",
    ylabel="Error Y [km]"
)
y=GM.lines!(axy, t.-tspan[1], errY)
axz = GM.Axis(
    fig[3,1],
    xlabel="Time[s]",
    ylabel="Error Z [km]"
)
y=GM.lines!(axz, t.-tspan[1], errZ)
fig