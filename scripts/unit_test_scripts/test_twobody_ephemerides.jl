using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

date = "2023-03-01T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-05T12:30:00"
endEpoch = utc2et(endDate)  # end date
a = twoBodyEarthEphemeride(epoch, endEpoch)

ν = get_heliocentric_position(a, endEpoch)
println(ν)

## Using the Ephemeride struct from Grant
tspan = (epoch, endEpoch)
eph = Ephemeride((tspan), 100, 399, 10, "J2000")

t = epoch + 64.9 # some random time\

state = interpolate(eph, t)

#unload kernels
unload("naif0012.tls")
unload("de440.bsp")