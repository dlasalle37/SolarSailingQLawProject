using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

###########################################################################################################################################################
# DEBUGGING:
# SETUP OF SCRIPT THAT CAUSED ERROR

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 20.0*86400 # amount of time [s] to simulate for
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
###########################################################################################################################################################
# VALUES THAT CAUSED ERROR
terr = 7.25849978740955e8;
oeerr = [
    9210.977469508176;
    0.19894763998216428;
    0.010001953085790594;
    -0.002032316953319921;
    89.99992476532344;
    1.6200477652642002
    ]
uerr =  [0.37695034132335264; 1.5700358587878527]
###########################################################################################################################################################

#RECREATE CONDITIONS
nu_earth = get_heliocentric_position(params.eph, terr)
params.current_time = terr
dxdt = [ # found in debugging 
    0.0072361175106722495;
    -3.511126653378224e-7;
    -1.9271342630806127e-10;
    4.1449959009927794e-6;
    2.018849499844806e-7;
    0.0007394410420640095
]
calculate_Q(oeerr, params) #Qcheck = 1.875537933485484e18