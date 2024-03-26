using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
import GLMakie as GM
using FrameTransformations

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
@time begin
    model = NormalizedGravityModel(n, m, l, R=6378.139, mu=398600.47);
    g = getFirstPartial(model, x, want_cf)
end
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
secondPartial = getSecondPartial(model, x, want_cf);

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
X0 = [
    -2750.0
    6848.0
    0.0
    -7.4732
    -2.9956
    1.0805
]

#Setup frame FrameTransformations
mdl = NormalizedGravityModel(n, m, l, R=6378, mu=398600.4418);
Orient.init_eop(datadir("iau2000a.eop.dat"))
FS = FrameSystem{4, Float64}()
@axes GCRF 1 GeocentricCelestialReferenceFrame
@axes ITRF 6 InternationalTerrestrialReferenceFrame
@point Earth 399
@point Spacecraft -1_900_000
@point acc -1_999_999
add_axes_inertial!(FS, GCRF)
add_point_root!(FS, Earth, GCRF)
add_axes_itrf!(FS, ITRF, GCRF) # adding the Earth-fixed reference
add_point_updatable!(FS, Spacecraft, Earth, GCRF) # adding the spacecraft as an updateable point in inertial coords
update_point!(FS, Spacecraft, X0, epoch)
add_point_updatable!(FS, acc, Earth, ITRF) # add a dummy point for the acceleration to rotate it later
update_point!(FS, acc, [0.0; 0.0; 0.0], epoch)


ps = (mu, FS, mdl)
prob = ODEProblem(two_body_eom_perturbed!, X0, (epoch, epoch+60000), ps, saveat=60)
sol = solve(prob, abstol=1.0e-12, reltol=1.0e-12)
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
    title = "Transfer A, Control computed within integration"
)

lin = GM.lines!(ax, orb[1,:], orb[2,:], orb[3,:], color=:limegreen, linewidth=2.0)

#write data
cart = hcat(sol.t, transpose(reduce(hcat, sol.u)))
open(datadir("ReportFile2.txt"),   "w") do io; writedlm(io,  cart); end