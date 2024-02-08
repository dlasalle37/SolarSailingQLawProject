using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))


open(datadir("kep.txt"), "r") do f
    global data = readdlm(f, '\t', Float64)
end

# plotting
xyz = Vector{Vector{Float64}}(undef, 0)
i = 1
for oe in eachrow(data)
    global mu = params.mu
    nue = get_heliocentric_position(eph, eph.t0+i*params.step_size)
    a = oe[1]; e = oe[2]; inc=oe[3]
    ω = oe[4]; RAAN = oe[5]+nue; θ = oe[6]
    (r, v) = coe2rv(a, e, inc, ω, RAAN, θ, mu)
    push!(xyz, r)
    global i+=1
end

x = []
y = []
z = []
for i in xyz
    push!(x, i[1])
    push!(y, i[2])
    push!(z, i[3])
end

plotlyjs()
plot(x, y, z)
xlabel!("Inertial X(km)"); ylabel!("Inertial Y (km)"); zlabel!("Inertial Z(km)")