using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

date = "2023-03-01T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-01T12:30:00"
endEpoch = utc2et(endDate)  # end date
a = twoBodyEarthEphemeride(epoch, endEpoch)

ν = get_heliocentric_position(a, endEpoch)
print(ν)

#unload kernels
unload("naif0012.tls")
unload("de440.bsp")