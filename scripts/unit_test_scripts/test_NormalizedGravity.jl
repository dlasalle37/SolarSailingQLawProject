using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM
using FrameTransformations
using BenchmarkTools

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP


l = datadir("EGM96_to360.ascii")
R = 6378
mu = 398600.4418
n = 5 # deg
m = 5 # order

x = [5489.150; 802.222; 3140.916] # appendix B vector from nasa technical doc
ep = x[3]/norm(x) # random value for epsilon, eps=z/r
want_cf = true

# ======== FIRST PARTIAL
Pg = normalized_legendre_generator(n,m, ep)
println("\nTiming $n by $m model creation")
model = @btime NormalizedGravityModel(n, m, l; R=6378.139, mu=398600.47);


println("\nTiming First Partial")
g = @btime getFirstPartial(model, x, want_cf)
println("\nAnalytic Acceleration [km/s2]:")
println(g)

# Calculating first partial via forwardDiff
#getPotential(model, x, want_cf)
dU(x) = getPotential(model, x, want_cf)
cfg = ForwardDiff.GradientConfig(dU, x) # Get the config
dUdx = ForwardDiff.gradient(dU, x, cfg)
diff = g - dUdx
println("\nDifference b/w analytical first partial and forwardDiff:")
println(diff)

# ======== SECOND PARTIAL
println("\nTiming Second Partial")
secondPartial = @btime getSecondPartial(model, x, want_cf);

# Calculating second partial (jacobian) with forwardDiff
d2U(x) = getFirstPartial(model, x, want_cf)
cfg = ForwardDiff.JacobianConfig(d2U, x) # Get the config
d2Udx2 = ForwardDiff.jacobian(d2U, x, cfg)
diff2 = secondPartial - d2Udx2
println("\nDifference b/w analytical second partial and forwardDiff:")
println(diff2)

# l2 = datadir("EGM2008_to2190_TideFree")
# mod2 = NormalizedGravityModelData(n, m, l2, R=6378.139, mu=398600.47, GravModel=:EGM2008);
# (g2, ~) = getFirstPartial(mod2, x, want_cf)

# ===== Running a simulation

# Setup Spacecraft
date = "2023-01-01T12:30:00" 
epoch = utc2et(date)
tof = 5*86400 # time of flight in seconds
coe0 = [10500.0, 0.150, 25*pi/180, 10.0*pi/180, 30.0*pi/180, 0.0] #[a, e, i, ape, ran, nu]
(r, v) = coe2rv(coe0[1], coe0[2], coe0[3], coe0[4], coe0[5], coe0[6], mu)
X0 = [r;v]

#Setup frame FrameTransformations
mdl = NormalizedGravityModel(n, m, l, R=6378.0, mu=398600.4418);
eop_load_data!(iers2010a, datadir("iau2000a.eop.dat"))
fs = FrameSystem{4, Float64}()
add_axes_icrf!(fs)
add_axes_gcrf!(fs)
add_axes_itrf!(fs, :ITRF, 23, 6)

ps = (mdl.mu, mdl, fs)
prob = ODEProblem(two_body_eom_perturbed!, X0, (epoch, epoch+tof), ps, saveat=60)
sol = solve(prob, Vern9(); abstol=1.0e-12, reltol=1.0e-12)
orb = reduce(hcat, sol.u)

# Plot 3d figure
fig = GM.Figure(;
)

ax = GM.Axis3(
    fig[1,1]; 
    aspect = :data, 
    xlabel = "x [km]", 
    ylabel = "y [km]", 
    zlabel = "z [km]",
    title = " "
)

lin = GM.lines!(ax, orb[1,:], orb[2,:], orb[3,:], color=:limegreen, linewidth=2.0)

#write data
cart = hcat(sol.t, transpose(reduce(hcat, sol.u)))
open(datadir("ReportFile2.txt"),   "w") do io; writedlm(io,  cart); end

fig
eop_unload_data!()