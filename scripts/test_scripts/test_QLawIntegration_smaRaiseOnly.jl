using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 120*86400 # amount of time [s] to simulate for
endTime = startTime+simTime

# QLaw Parameter setup
eph = twoBodyEarthEphemeride(startTime, endTime)  # create the earth ephemeride
sc = basicSolarSail()
X0 = [9222.7; 0.20; 0.573*pi/180; 0.00; 90; 0.0]  # COE initial conditions [a, e, i, argPer, RAAN, trueAnom]
XT = [13000.0, 0.75, 0.01*pi/180, 270.0*pi/180, 90*pi/180] # Targets # note that targets has 5 elements, while X0 has 6
oetols = [100, 0.001, 0.01, 0.01, 0.01]
Woe = [1.0, 0.0, 0.0, 0.0, 0.0]
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
    writeData=true,
    kimp=100
    )

#######################################################################################################################################################
#Set up and solve ode
abstol = 1.0E-6
reltol = 1.0E-6
tspan = (startTime, endTime)
prob = ODEProblem(QLawEOM!, X0, tspan, params, abstol=abstol, reltol=reltol)
condition(u, t, integrator) = maximum([
    params.Woe[1]*abs(u[1] - params.oet[1]) - params.oeTols[1];
    params.Woe[2]*abs(u[2] - params.oet[2]) - params.oeTols[2];
    params.Woe[3]*abs(u[3] - params.oet[3]) - params.oeTols[3];
    params.Woe[4]*abs(acos(cos(u[4] - params.oet[4]))) - params.oeTols[4];
    params.Woe[5]*abs(acos(cos(u[5] - params.oet[5]))) - params.oeTols[5]
]) <= 0
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)
sol = solve(prob,  callback=cb, saveat=60)

print("End Values: ")
println(sol.u[end])
#######################################################################################################################################################
# Writing
if params.writeData
    open(datadir("kep.txt"),   "w") do io; writedlm(io,  sol.u); end
end

# ===== Plotting
# First read the data
kep = readdlm(datadir("kep.txt"), '\t', '\n'; header=false)

# Convert to cartesian
cart = Matrix{Float64}(undef, size(kep))
for row in axes(kep, 1)
    r, v = coe2rv(kep[row,1], kep[row,2], kep[row,3], kep[row,4], kep[row,5], kep[row,6], 398600.4418)
    cart[row,1:3] .= r
    cart[row,4:6] .= v
end

# Plot 3d figure
fig = GM.Figure(;
    size = (1920,1080),
)

ax = GM.Axis3(
    fig[1,1]; 
    aspect = :data, 
    xlabel = "x [km]", 
    ylabel = "y [km]", 
    zlabel = "z [km]",
)

GM.lines!(ax, cart[:,1], cart[:,2], cart[:,3], color=:blue, linewidth=0.5)

# Unload Kernels
unload("naif0012.tls")
unload("de440.bsp")