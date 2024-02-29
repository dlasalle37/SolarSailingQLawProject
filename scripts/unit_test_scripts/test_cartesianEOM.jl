using DrWatson

@quickactivate "SolarSailingQLawProject"

# Here you may include files from the source directory
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP

sc=basicSolarSail();

# Set up condition
date = "2023-03-01T12:30:00" 
startTime = utc2et(date)  # start date
simTime = 30*24*3600.0
endTime = startTime + simTime
tspan = (0, simTime)
eph = twoBodyEarthEphemeride(startTime, endTime)
p = (398600.4418, sc, eph);  # parameters tuple as ordered in the EOM
X0 = [42164.0;0.0;0.0;0.0;4.03;0.00]  # initial position and velocity vector

# Tolerances
reltol=1.0E-6
abstol=1.0E-6

prob = ODEProblem(solarSailEOM_cartesian!, X0, tspan, p, abstol=abstol, reltol=reltol)
sol = solve(prob);
plotlyjs()
plot(sol, idxs=(1,2,3))
title!("Start  Date: $date")
xlabel!("Inertial X(km)"); ylabel!("Inertial Y (km)"); zlabel!("Inertial Z(km)")