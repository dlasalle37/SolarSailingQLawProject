using DrWatson
using SPICE
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

mu = 398600.4418
r1 = [42164.0, 36.4, 0.971] # km
v1 = [0.21, 3.031, 0.19] # km/sail

coe = rv2coe(r1, v1, mu)

r2, v2 = coe2rv(coe[1], coe[2], coe[3], coe[4], coe[5], coe[6], mu)

println("initial r: $r1 \ninitial v: $v1")
println("check r: $r2 \ncheck v: $v2")