using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("utils.jl"))
using BenchmarkTools
mu = 398600.4418
r = [40682069.1092143, -11057054.204445519, 4373.164482962849]*1e-3
v = [806.5489379693863, 2967.530987081164, -1.7511311698939385]*1e-3

rr = [40695619.3157, -10996325.2888, 16862.675]*1e-3

mat = @btime eci2ric(r, v)
ans=mat*(rr-r)

expect = [-2.85066, 62.14974, 12.525199]
println(ans-expect)