using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Simulation time setup:
date = "2023-03-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 150.0*24*3600 # amount of time [s] to simulate for
endTime = startTime+simTime

eph = twoBodyEarthEphemeride(startTime, endTime)  # create the earth ephemeride
trueAnom_earth_i = eph.trueAnom_initial
lambda_i = 0.0-trueAnom_earth_i
sc = basicSolarSail(C1=2.0 ,C2=0.0, C3=0.0, areaParam=2.5E-5)
mu = 398600.4418;
X0 = [42164.0; 0.00; 0.001*pi/180; 0.00; 0.00; 0.0]  # COE initial conditions [a, e, i, argPer, RAAN, trueAnom]
XT = [42764.0, 0.7, 20*pi/180, 270.0*pi/180, 90*pi/180] # Targets # note that targets has 5 elements, while X0 has 6
oetols = [10, 0.01, 0.1, 0.1, 0.1]
p = QLawParams(sc, eph, X0, XT, oetols)
p.alpha = 30*pi/180
p.beta = 0.0*pi/180
tspan = (0, simTime)

# Tolerances
reltol=1.0E-6
abstol=1.0E-6
prob = ODEProblem(gauss_variational_eqn!, X0, tspan, p, abstol=abstol, reltol=reltol, saveat=600)
sol = solve(prob);

# Convert all to r, v:
solCartesian = zeros(length(sol.t), 6)
k=1
for i in sol.u
    trueAnom_earth = get_heliocentric_position(eph, startTime+sol.t[k])
    a = i[1]; e = i[2]; inc=i[3]
    ω = i[4]; RAAN = i[5]+trueAnom_earth; θ = i[6]
    (r, v) = coe2rv(a, e, inc, ω, RAAN, θ, mu)
    solCartesian[k,:] = [r;v]
    global k+=1
end

plotlyjs()
plot(solCartesian[:,1], solCartesian[:,2], solCartesian[:,3])
xlabel!("Inertial X(km)"); ylabel!("Inertial Y (km)"); zlabel!("Inertial Z(km)")