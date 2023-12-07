###################################################################################################################################################
# Set up test conditions
using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

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

a = 9210.977469508176
e = 0.19894763998216428
inc = 0.010001953085790594
ape = -0.002032316953319921
lam =  89.99992476532344
tru = 1.6200477652642002
nustar_a = 2.6790707975069004
sig_a = -1.0
eps_a = [1; 0; 0; 0; 0; 0]
tru_E = 6.21595532563805
f0 = [0.0; 0.0; 0.0; 0.0; -2.0586466118771686e-7; 0.0007439937671808995]
###################################################################################################################################################
# Call function
adotnn_test = oedotnn(a, e, inc, ape, lam, tru, nustar_a, sig_a, eps_a, tru_E, f0, params)
println("Function Result: adotnn = $adotnn_test")