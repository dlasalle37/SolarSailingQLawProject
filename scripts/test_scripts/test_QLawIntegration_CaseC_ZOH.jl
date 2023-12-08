using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
using DelimitedFiles
import GeometryBasics as GB
import GLMakie as GM

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 45.0*24*3600 # amount of time [s] to simulate for
endTime = startTime+simTime

# QLaw Parameter setup
eph = twoBodyEarthEphemeride(startTime, endTime)  # create the earth ephemeride
sc = basicSolarSail()
X0 = [9222.7; 0.20; 0.573*pi/180; 0.00; 90; 0.0]  # COE initial conditions [a, e, i, argPer, RAAN, trueAnom]
XT = [26500.0, 0.75, 0.01*pi/180, 270.0*pi/180, 90*pi/180] # Targets # note that targets has 5 elements, while X0 has 6
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
    writeData=true,
    kimp=100
    )
finalOE, exitcode = QLawIntegrator(params)

println("Starting Values: ", X0[1:6])
println("Final Values: ", finalOE[1:6])
println("Targets: ", XT)
println(exitcode)

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

#Pull starting, ending points
startPoint = cart[1, 1:3]
endPoint = cart[end, 1:3]

# Plot 3d figure
fig = GM.Figure(;
    size = (1920,1080)
)

ax = GM.Axis3(
    fig[1,1]; 
    aspect = :data, 
    xlabel = "x [km]", 
    ylabel = "y [km]", 
    zlabel = "z [km]",
    title = "Apoapsis-Raising Transfer computed with QLaw"
)

lin = GM.lines!(ax, cart[:,1], cart[:,2], cart[:,3], color=:blue, linewidth=0.5)

# Plot start/end points
sP = GM.scatter!(ax, startPoint[1], startPoint[2], startPoint[3], markersize=10.0, color=:black)
eP = GM.scatter!(ax, endPoint[1], endPoint[2], endPoint[3], markersize=10.0, color=:red)

# Create and add a sphere to represent earth
sphere = Sphere(Point3f(0), 6378.0)
spheremesh = GB.mesh(Tesselation(sphere, 64))
sph = GM.mesh!(ax, spheremesh; color=(:blue))

#Create legend
GM.Legend(fig[1, 2], [lin, sP, eP], ["Satellite Trajectory", "Starting Point", "Ending Point"])

    # Unload kernels
unload("naif0012.tls")
unload("de440.bsp")