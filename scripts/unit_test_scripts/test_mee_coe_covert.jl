using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

# Check conversions
kep = [6378.0, 0.15, pi/4, pi/6, pi/6, pi/2]
(p, f, g, h, k, L) = coe_to_mee(kep...)
(a, e, i, Ω, ω, θ) = mee_to_coe(p, f, g, h, k, L)
diff = kep - [a, e, i, Ω, ω, θ]
println("Difference cause by converting back and forth (should be zeros): $diff")

# Check MEE vs COE dynamics