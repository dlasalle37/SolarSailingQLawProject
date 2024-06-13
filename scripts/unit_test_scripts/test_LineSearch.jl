using DrWatson
@quickactivate("SolarSailingQLawProject")
include(srcdir("LineSearch.jl"))
import GLMakie as gm

fn(x) = exp(0.1*x)*sin(x)
mn = -5
mx = 5

xstar, fxstar = gss(fn, mn, mx)

xs = collect(mn:.1:mx)
ys = Vector{Float64}(undef, 0)
for x in xs
    push!(ys, fn(x))
end

fig = gm.Figure()
ax = gm.Axis(fig[1,1])
gm.lines!(ax, xs, ys)
gm.scatter!(ax, xstar, fxstar, color=:Red)
fig