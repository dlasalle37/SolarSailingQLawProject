using DrWatson
using SPICE
@quickactivate "SolarSailingQLawProject"
include(srcdir("utils.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

date = "2023-06-01T012:00:00" 
et = utc2et(date)  # start date

tspan = 0:600:360.0*24*3600

k=1
r = zeros(length(collect(tspan)), 3)
for t in tspan
    (X, o) = spkez(399, et+t, "ECLIPJ2000", "none", 10);
    rES_vec = X[1:3]; vES_vec = X[4:6]  # position and velocity of Earth wrt Sun in J2000 Ecliptic frame
    #earth_coe = rv2coe(rES_vec, vES_vec, 1.327E11)  # earth's keplerian orbital elements in its heliocentric orbit
    #(r[i,:], ~) = coe2rv
    r[k, :] = rES_vec
    global k+=1
end



plot(r[:,1], r[:,2], r[:,3])
xlabel!("Inertial X(km)"); ylabel!("Inertial Y (km)"); zlabel!("Inertial Z(km)")