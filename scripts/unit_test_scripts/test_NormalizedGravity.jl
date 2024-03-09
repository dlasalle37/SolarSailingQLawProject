using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

l = datadir("EGM96_to360.ascii")
R = 6378
mu = 398600
n = 5 # deg
m = 5 # order

x = [5489.150; 802.222; 3140.916] # appendix B vector from nasa technical doc
ep = x[3]/norm(x) # random value for epsilon, eps=z/r
want_cf = true

@time begin
    model = NormalizedGravityModelData(n, m, l, R=6378.139, mu=398600.47);
    (g, P) = getFirstPartial(model, x, want_cf)
    Pg = normalized_legendre_generator(n, m, ep); # checking legendre calcs
    println(P==Pg)
end
println("Acceleration [km/s]:")
println(g)

d2U = getSecondPartial(model, x, want_cf);

l2 = datadir("EGM2008_to2190_TideFree")
mod2 = NormalizedGravityModelData(n, m, l2, R=6378.139, mu=398600.47, GravModel=:EGM2008);
(g2, ~) = getFirstPartial(mod2, x, want_cf)