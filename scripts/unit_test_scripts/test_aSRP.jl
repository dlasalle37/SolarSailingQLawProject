# Setup
using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 20*86400 # amount of time [s] to simulate for
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

# arbitrary control params
u = [pi/4; pi/8]
nue = eph.trueAnom_initial
aSRP_test = aSRP(u, sc, eph, nue)


# Check 
a = u[1]
b = u[2]
C1 = sc.C[1]
C2 = sc.C[2]
C3 = sc.C[3]
a_over_m = sc.areaParam
d = distance_to_sun(eph, nue)
G0 = 1.02E14

aSRP_check = a_over_m*G0/d^2 * cos(a)*[
    C1*cos(a)^2+C2*cos(a) + C3;
    -((C1*cos(a)+C2)*sin(a)*sin(b));
    -((C1*cos(a)+C2)*sin(a)*cos(b))
]

println(aSRP_check)
println(aSRP_test)
println("difference: $(aSRP_check-aSRP_test)")