using DrWatson
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

## Using the Ephemeride struct from Grant
tspan = (epoch, endEpoch)
eph = Ephemeride((tspan), 1000, 399, 10, "ECLIPJ2000") # targ is earth, obs is sun, so vector points from sun->earth 

t = epoch + 64.4 # some random time\

state = getState(eph, t)

eph2 = Ephemeride((tspan), 1000, 10, 399, "ECLIPJ2000") # targ is sun
state2 = getState(eph2, t)

println(state)
println(state2)

# Testing the eclipse calculator
scpos = 8000 * [-0.9416839528031331
0.3364986371609293
-1.4896439549352508e-5]
println(isEclipsed(eph, scpos, t))


#unload kernels
unload("naif0012.tls")
unload("de440.bsp")