using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
include(srcdir("ThirdBodyPerturbations.jl"))
import GLMakie as GM
import GeometryBasics as GB

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP

# Simulation time setup:
mu = 398600.4418
date = "2023-01-01T12:30:00" 
startTime = utc2et(date)  # start date in seconds past j2000
simTime = 10*86400 # amount of time [s] to simulate for
endTime = startTime+simTime

# QLaw Parameter setup
eph = Ephemeride((startTime, endTime), 1000, 399, 10, "J2000")

# Third body perturbations info
perturbing_body = eph.obs # taking the sun as the observer

X0 = [42164.0, 0.0, 0.0, 0.0, -3.03, 0.0]
ps = (mu, eph, perturbing_body)
abstol = 1.0E-6
reltol = 1.0E-6
tspan = (startTime, endTime)
prob = ODEProblem(two_body_eom_perturbed_tb!, X0, tspan, ps, abstol=abstol, reltol=reltol)
sol = solve(prob, Vern9(), saveat=60)
t = sol.t
cart = reduce(hcat, sol.u)';
kep = Matrix{Float64}(undef, size(cart))
for row in axes(cart, 1)
    r = cart[row, 1:3]
    v = cart[row, 4:6]
    coe = rv2coe(r, v, mu)
    kep[row,:] = [coe[1], coe[2], coe[3], coe[4], coe[5], coe[6]]
end

# Plot 3d figure
fig = GM.Figure(;
)
ax = GM.Axis3(
    fig[1,1]; 
    aspect = :data, 
    xlabel = "x [km]", 
    ylabel = "y [km]", 
    zlabel = "z [km]",
    title = "Third body perturbations"
)
lin = GM.lines!(ax, cart[:,1], cart[:,2], cart[:,3], color=:blue, linewidth=0.5)
## Create and add a sphere to represent earth
sphere = GB.Sphere(GB.Point3f(0), 6378.0)
spheremesh = GB.mesh(GB.Tesselation(sphere, 64))
sph = GM.mesh!(ax, spheremesh; color=(:blue))

#plot oe histories
fig2 = GM.Figure(title="Orbital Element Histories")
axa = GM.Axis(
    fig2[1,1],
    xlabel="Time[days]",
    ylabel="km",
    title="Semi-Major Axis",
)
axe = GM.Axis(
    fig2[2,1],
    xlabel="Time[days]",
    title="Eccentricity",
)
axi = GM.Axis(
    fig2[1,2],
    xlabel="Time[days]",
    ylabel="Angle[Deg]",
    title="Inclination",
)
axape = GM.Axis(
    fig2[2,2],
    xlabel="Time[days]",
    ylabel="Angle[Deg]",
    title="Arg.Perigee",
)
axRAN  = GM.Axis(
    fig2[3,1],
    xlabel="Time[days]",
    ylabel="Angle[deg]",
    title="RAAN",
)
axtru = GM.Axis(
    fig2[3,2],
    xlabel="Time[days]",
    ylabel="Angle[deg]",
    title="True anom",
)


GM.lines!(axa, t/86400, kep[:,1])
GM.lines!(axe, t/86400, kep[:,2])
GM.lines!(axi, t/86400, kep[:,3]*180/pi)
GM.lines!(axape, t/86400, kep[:,4]*180/pi)
GM.lines!(axRAN, t/86400, kep[:,5]*180/pi)
GM.lines!(axtru, t/86400, kep[:,6]*180/pi)