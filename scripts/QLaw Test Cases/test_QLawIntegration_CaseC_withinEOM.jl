using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM
import GeometryBasics as GB

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP

# Simulation time setup:
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 50*86400 # amount of time [s] to simulate for
endTime = startTime+simTime

# QLaw Parameter setup
eph = Ephemeride((startTime, endTime), 1000, 399, 10, "J2000")
sc = basicSolarSail()
coee = getCOE(eph, eph.t0)
nue = coee[6]
X0 = [26500; 0.70; 0.573*pi/180; 0.000; 2.0354056994857928-nue; 0.0]  # COE initial conditions [a, e, i, argPer, lambda, trueAnom]
XT = [26700, 0.75, 0.2*pi/180, 30*pi/180, 90.0*pi/180] # Targets # note that targets has 5 elements, while X0 has 6
oetols = [10, 0.001, 0.01, 0.01, 0.01]
Woe = [1.0, 1.0, 1.0, 1.0, 0.0]
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
@time begin 
    sol = solve(prob, AutoTsit5(Rosenbrock23()),  saveat=60, callback=ccb) 
end
print("End Values: ")
println(sol.u[end])
#######################################################################################################################################################
# ===== Post-Processing
x = reduce(hcat, sol.u)' # full solution (matrix form)
t = sol.t.-params.eph.t0 # shift time back to start at zero

kep = Matrix{Float64}(undef, size(x))
cart = Matrix{Float64}(undef, size(x))
lambda = Vector{Float64}(undef, size(x)[1])
angles = Matrix{Float64}(undef, (size(x)[1], 2))
shadex = []
shadey = [] 
shadez = []
for row in axes(x, 1)
    local coe = getCOE(eph, eph.t0+t[row])
    local nue = coe[6] # pull true anomaly

    # Save lambda separately
    lambda[row] = x[row, 5]

    # Convert to keplerian orbital elements (from lambda = RAAN-nue)
    kep[row,:] = [x[row,1], x[row,2], x[row,3], x[row,4], x[row,5]+nue, x[row,6]]

    # Convert keplerian oe's to cartesian coords
    r, v = coe2rv(x[row,1], x[row,2], x[row,3], x[row,4], x[row,5]+nue, x[row,6], 398600.4418)
    cart[row,1:3] .= r
    cart[row,4:6] .= v

    # Gather Eclipse info
    local eclipsed = isEclipsed(eph, r, eph.t0+t[row])
    if eclipsed
        push!(shadex, cart[row, 1])
        push!(shadey, cart[row, 2])
        push!(shadez, cart[row, 3])
    end


    # Re-compute the sail angles
    local alpha, beta, dQdx = compute_control(kep[row,:], params)
    angles[row,:] = [alpha, beta]
end

# Writing PP Output
if params.writeData
    # Writing to the datadir, paste these files into scripts\data for matlab plotting
    open(datadir("kep.txt"),   "w") do io; writedlm(io,  kep); end # keplerian oe set
    open(datadir("lambda.txt"), "w") do io; writedlm(io, lambda); end # oe lambda
    open(datadir("eclipsed.txt"), "w") do io; writedlm(io, [shadex shadey shadez]); end # the eclipsed points
    open(datadir("cart.txt"), "w") do io; writedlm(io, cart); end # full cartesian trajectory
    open(datadir("angles.txt"), "w") do io; writedlm(io, angles); end # sail angles
    open(datadir("discrete_times.txt"), "w") do io; writedlm(io, t); end # time vector to plot against
end



# ===== Plotting
# First read the data
kep = readdlm(datadir("kep.txt"), '\t', '\n'; header=false)
cart = readdlm(datadir("cart.txt"), '\t', '\n'; header=false)

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
sP = GM.scatter!(ax, startPoint[1], startPoint[2], startPoint[3], markersize=10.0, color=:green)
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
GM.lines!(axlam, t/86400, lambda*180/pi)
GM.lines!(axRAN, t/86400, kep[:,5]*180/pi)
GM.lines!(axtru, t/86400, kep[:,6]*180/pi)




# Unload Kernels
unload("naif0012.tls")
unload("de440.bsp")