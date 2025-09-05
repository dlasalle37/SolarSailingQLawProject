using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
using BenchmarkTools

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP

sc = basicSolarSail()
X0 = [26500; 0.70; 0.573*pi/180; 0.000; -4.179809375168203; 0.0]
XT = [26700, 0.75, 0.2*pi/180, 30*pi/180, 90.0*pi/180] 
oetols = [10, 0.001, 0.01, 0.01, 0.01]
Woe = [1.0, 1.0, 1.0, 1.0, 0.0]
qlaw_type = Oguri 

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 155*86400 # amount of time [s] to simulate for
endTime = startTime+simTime

# Gravity model
n = 5
m = 5
l = datadir("EGM96_to360.ascii")
mdl = NormalizedGravityModel(n, m, l, R=6378.139, mu=398600.4418);

# ======= Frame System Setup
eop_load_data!(iers2010a, datadir("iau2000a.eop.dat"))
fs = FrameSystem{4, Float64}()
add_axes_icrf!(fs)
add_axes_gcrf!(fs)
add_axes_itrf!(fs, :ITRF, 23, 6)

# ======= QLaw Parameter setup
eph = Ephemeride((startTime, endTime), 1000, 399, 10, "J2000")
params = QLawParams(
    sc,
    eph, 
    X0, 
    XT, 
    oetols,
    fs,
    mdl,
    Woe=Woe,
    rp_min=6578.0,
    a_esc=1.0E5,
    max_sim_time = simTime,
    step_size = 60.0,
    kimp=100,
    writeData=true,
    type=qlaw_type
    )

# Calculate partials
X = X0
t = eph.t0

Q(x) = calculate_Q(t, x, params) # defining a unary closure to allow for passage of params into ForwardDiff
cfg = ForwardDiff.GradientConfig(Q, X) # Get the config

println("AD timing:")
@btime dQdx_fd = ForwardDiff.gradient(Q, X, cfg)

println("Mixed analytical/AD timing")
@btime dQdx_mixed = calculate_semi_analytic_partials(t, X, params)

dQdx_fd = ForwardDiff.gradient(Q, X, cfg)
dQdx_mixed = calculate_semi_analytic_partials(t, X, params)
println(dQdx_fd-dQdx_mixed)

# Unload Kernels
unload("naif0012.tls")
unload("de440.bsp")
eop_unload_data!()