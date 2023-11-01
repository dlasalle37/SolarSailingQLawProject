using DrWatson
@quickactivate "SolarSailingQLawProject"

include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

epoch = "2018-06-06T20:45:00"
eti = utc2et(epoch)
simTime =3*360*24*3600
etf = eti+simTime
tspan = 1:24*3600:simTime
ephemeride = twoBodyEarthEphemeride(eti, etf)
k=1
dist = zeros(length(collect(tspan)))
for i in tspan
    nu = get_heliocentric_position(ephemeride, eti+i)
    dist[k] = distance_to_sun(ephemeride, nu)
    global k+=1
end

plot(collect(tspan)./(24*3600), dist)
xlabel!("Time(days)"); ylabel!("Distance(km)")
title!("Distance to sun (Epoch: $epoch)")
