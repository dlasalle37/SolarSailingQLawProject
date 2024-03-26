using DrWatson
using BenchmarkTools
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP


date = "2023-03-01T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-05T12:30:00"
endEpoch = utc2et(endDate)  # end date

## Using the Ephemeride
t = epoch + 64.4 # some random time\
tspan = (epoch, endEpoch)
eph = Ephemeride((tspan), 1000, 399, 10, "ECLIPJ2000") # targ is earth, obs is sun, so vector points from sun->earth 

scpos = 8000 * [-0.9416839528031331
0.3364986371609293
-1.4896439549352508e-5]
@benchmark isEclipsed(eph, scpos, t) setup=(eph=eph, scpos=scpos, t=t)