using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM
import GeometryBasics as GB

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 55*86400 # amount of time [s] to simulate for
endTime = startTime+simTime

# ======= QLaw Parameter setup
eph = twoBodyEarthEphemeride(startTime, endTime)  # create the earth ephemeride
sc = basicSolarSail()
nue = get_heliocentric_position(eph, eph.t0)
X0 = [9222.7; 0.20; 0.573*pi/180; 0.00; 2.02; 0.0]  # COE initial conditions [a, e, i, argPer, lamba, trueAnom]
XT = [26500.0, 0.75, 0.01*pi/180, 0.0*pi/180, 90*pi/180] # Targets # note that targets has 5 elements, while X0 has 6
oetols = [10, 0.001, 0.01, 0.01, 0.01]
Woe = [1.0, 1.0, 1.0, 0.0, 0.0]
params = QLawParams(
    sc, 
    eph, 
    X0, 
    XT, 
    oetols, 
    Woe=Woe,
    rp_min=6578.0,
    a_esc=1.0E5,
    max_sim_time = simTime,
    step_size = 60.0,
    kimp=100,
    writeData=true
    )

# ======= Integrator Setup
abstol = 1.0E-6
reltol = 1.0E-6
tspan = (startTime, endTime)
prob = ODEProblem(QLawEOM!, X0, tspan, params, abstol=abstol, reltol=reltol)
# ======= Callback Setups for Integrator:
# ======= To use a Discrete callback, comment in this block
#=condition(u, t, integrator) = callback_function_error_check(u, t, params)
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)=#

# ======= To use a Continuous callback, comment in this block
condition(u, t, integrator) = continuous_callback_errorcheck(u, t, params)
affect!(integrator) = terminate!(integrator)
ccb = ContinuousCallback(condition, affect!)

# ====== Run solve function to solve DE
sol = solve(prob, AutoTsit5(Rosenbrock23()),  saveat=600, callback=ccb)

print("End Values: ")
println(sol.u[end])
#######################################################################################################################################################
# Writing
if params.writeData
    open(datadir("kep.txt"),   "w") do io; writedlm(io,  sol.u); end
end

# ===== Post-Processing and Plotting
# First read the data
kep = readdlm(datadir("kep.txt"), '\t', '\n'; header=false)
t = sol.t.-params.eph.t0 # shift time back to start at zero

# Some Post-processing items
cart = Matrix{Float64}(undef, size(kep))
ran = Vector{Float64}(undef, size(kep)[1])
angles = Matrix{Float64}(undef, (size(kep)[1], 2))
for row in axes(kep, 1)
    local nue = get_heliocentric_position(eph, eph.t0+t[row]) # get earth pos @ each time for conversion to keplerian elemetns
    ran[row] = kep[row,5]+nue # save RAAN to its own variable to plot later

    # Get cartesian Coords
    r, v = coe2rv(kep[row,1], kep[row,2], kep[row,3], kep[row,4], kep[row,5]+nue, kep[row,6], 398600.4418)
    cart[row,1:3] .= r
    cart[row,4:6] .= v

    # Re-compute the sail angles
    alpha, beta, dQdx = compute_control(kep[row,:], params)
    angles[row,:] = [alpha, beta]
end

#Pull starting, ending points
startPoint = cart[1, 1:3]
endPoint = cart[end, 1:3]

# Plot 3d figure
fig = GM.Figure(;
)

ax = GM.Axis3(
    fig[1,1]; 
    aspect = :data, 
    xlabel = "x [km]", 
    ylabel = "y [km]", 
    zlabel = "z [km]",
    title = "Transfer A, Control computed within integration"
)

lin = GM.lines!(ax, cart[:,1], cart[:,2], cart[:,3], color=:blue, linewidth=0.5)

# Plot start/end points
sP = GM.scatter!(ax, startPoint[1], startPoint[2], startPoint[3], markersize=10.0, color=:black)
eP = GM.scatter!(ax, endPoint[1], endPoint[2], endPoint[3], markersize=10.0, color=:red)

## Create and plot initial/final orbits 
#Initial:
mu = params.mu
X0 = cart[1,:]
a0 = kep[1,1]
period_initial = 2*pi/sqrt(mu/a0^3)
prob = ODEProblem(two_body_eom!, X0, (0, period_initial), mu, saveat=60)
sol2 = solve(prob)
orb0 = reduce(hcat, sol2.u)
lin2 = GM.lines!(ax, orb0[1,:], orb0[2,:], orb0[3,:], color=:limegreen, linewidth=2.0)
        
#Final:
af = kep[end, 1]
XF = cart[end, :]
period_final = 2*pi/sqrt(mu/af^3)
prob = ODEProblem(two_body_eom!, XF, (0, period_final), mu, saveat=60)
sol3 = solve(prob)
orbF = reduce(hcat, sol3.u)
lin3 = GM.lines!(ax, orbF[1,:], orbF[2,:], orbF[3,:], color=:red, linewidth=2.0)
    
    
## Create and add a sphere to represent earth
sphere = GB.Sphere(GB.Point3f(0), 6378.0)
spheremesh = GB.mesh(GB.Tesselation(sphere, 64))
sph = GM.mesh!(ax, spheremesh; color=(:blue))
    
#Create legend
GM.Legend(fig[1, 2], [lin, sP, eP, lin2, lin3], ["Satellite Trajectory", "Starting Point", "Ending Point", "Initial Orbit", "Final Orbit"])

# Plot steering Law
fig2 = GM.Figure()
ax2 = GM.Axis(
    fig2[1,1];
    xlabel = "Time[days]",
    ylabel = "Angle [Deg]",
    title = "Steering Angle [alpha]"
    )
alp = GM.lines!(ax2, t/86400, angles[:,1]*180/pi)
ax3 = GM.Axis(
    fig2[2,1];
    xlabel = "Time[days]",
    ylabel = "Angle [Deg]",
    title = "Steering Angle [beta]"
    )
bet = GM.lines!(ax3, t/86400, angles[:,2]*180/pi)

#plot oe histories
fig3 = GM.Figure(title="Orbital Element Histories")
axa = GM.Axis(
    fig3[1,1],
    xlabel="Time[days]",
    ylabel="km",
    title="Semi-Major Axis",
)
axe = GM.Axis(
    fig3[2,1],
    xlabel="Time[days]",
    title="Eccentricity",
)
axi = GM.Axis(
    fig3[1,2],
    xlabel="Time[days]",
    ylabel="Angle[Deg]",
    title="Inclination",
)
axape = GM.Axis(
    fig3[2,2],
    xlabel="Time[days]",
    ylabel="Angle[Deg]",
    title="Arg.Perigee",
)
axlam = GM.Axis(
    fig3[3,1],
    xlabel="Time[days]",
    ylabel="Angle[Deg]",
    title="Lambda",
)
axtru = GM.Axis(
    fig3[3,3],
    xlabel="Time[days]",
    ylabel="Angle[deg]",
    title="True anom",
)
axRAN  = GM.Axis(
    fig3[3,2],
    xlabel="Time[days]",
    ylabel="Angle[deg]",
    title="RAAN",
)

GM.lines!(axa, t/86400, kep[:,1])
GM.lines!(axe, t/86400, kep[:,2])
GM.lines!(axi, t/86400, kep[:,3]*180/pi)
GM.lines!(axape, t/86400, kep[:,4]*180/pi)
GM.lines!(axlam, t/86400, kep[:,5]*180/pi)
GM.lines!(axRAN, t/86400, ran*180/pi)
GM.lines!(axtru, t/86400, kep[:,6]*180/pi)


# Unload Kernels
unload("naif0012.tls")
unload("de440.bsp")