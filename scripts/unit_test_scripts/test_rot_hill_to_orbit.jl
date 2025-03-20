# Set up test conditions
using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

x =  [9210.977469508176;
0.19894763998216428;
0.010001953085790594;
-0.002032316953319921;
89.99992476532344;
1.6200477652642002;] # radians

inc = x[3]; ape = x[4]; lam = x[5]; tru = x[6];
res_test = hill_to_orbit_transform(i, ape, lam, tru)

# Check
ang = tru+ape
Mat1 = [
    cos(ang) sin(ang) 0;
    -sin(ang) cos(ang) 0;
    0 0 1
]
Mat2 = [
    1 0 0;
    0 cos(inc) sin(inc);
    0 -sin(inc) cos(inc)
]
Mat3 = [
    cos(lam) sin(lam) 0;
    -sin(lam) cos(lam) 0;
    0 0 1
]
res_check = Mat1*Mat2*Mat3

println(res_test)
println(res_check)
res_test == res_check