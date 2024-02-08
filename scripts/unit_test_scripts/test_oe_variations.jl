############################################################################################################################################
# Set up test conditions
using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

nudot_earth = 1.0E-7
a = 15000.0 #km
e = 0.15
inc = pi/6;
ape = pi/4
lam = 0.0
tru = 30*pi/180
mu = 398600.4418
##########################################################################################################################################

# Call the two functions that calculate variations 
# (note QLaw has its own internal function F to calculate the matrix to keep everything encapsulated)
tol = 1.0E-9
F1 = F(a, e, inc, ape, tru, mu) # QLaw's Calculate
(f0, F2) = augmented_keplerian_varaitions(a, e, inc, ape, lam, tru, nudot_earth, mu)
check1 = F2-F1