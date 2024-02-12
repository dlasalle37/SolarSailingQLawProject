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

## Using the Ephemeride struct from Grant
tspan = (epoch, endEpoch)
eph = Ephemeride((tspan), 1000, 399, 10, "ECLIPJ2000")

t = epoch + 64.4 # some random time\

state = getState(eph, t)

#unload kernels
unload("naif0012.tls")
unload("de440.bsp")