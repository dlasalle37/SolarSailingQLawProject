using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

l = datadir("testmat.txt")
R = 6378
mu = 398600
nm = 300

#open(datadir("testmat.txt"), "w") do io; writedlm(io, rand(nm, nm)); end # create random nxm matrix
x = [3200; -500; 6378]
want_cf = true


model = NormalizedGravityModelData(nm, l, l);
(g, P) = getFirstPartial(model, x, want_cf)


ep = 0.8916235738552135 # random value for epsilon, eps=z/r
Pg = normalized_legendre_generator(nm, nm, ep); # checking legendre calcs
P==Pg