using DrWatson
@quickactivate("SolarSailingQLawProject")
include(srcdir("LineSearch.jl"))
import GLMakie as gm
import Optim as opt
using BenchmarkTools

fn(x) = exp(0.1*x)*sin(x)
mn = -5
mx = 5

xstar, fxstar, iter = gss(fn, mn, mx)
println("Custom golden search converged in $iter iterations")
xs = collect(mn:.1:mx)
ys = Vector{Float64}(undef, 0)
for x in xs
    push!(ys, fn(x))
end
 
sol_opt = opt.optimize(fn, mn, mx, opt.GoldenSection(), eps=1e-5)
iter_opt = sol_opt.iterations
minimizer_opt = sol_opt.minimizer
min_opt = sol_opt.minimum
println("Optim.jl Golden search converged in $iter_opt iterations")

println("\nTIMING COMPARISON")
println("\nCustom method:")
@btime gss(fn, mn, mx)

println("\nOptim.jl:")
@btime opt.optimize(fn, mn, mx, opt.GoldenSection(), eps=1e-5)

fig = gm.Figure()
ax = gm.Axis(fig[1,1])
gm.lines!(ax, xs, ys, label="Function")
gm.scatter!(ax, xstar, fxstar, color=:Red, markersize=15, label="Custom GSS")
gm.scatter!(ax, minimizer_opt, min_opt, marker=:cross, color=:Blue, label="Optim.jl GSS")
gm.axislegend(ax)
fig

# println("\nBenchmarking Custom GSS")
# @benchmark gss(fn, mn, mx)