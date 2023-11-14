using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Spacecraft/epoch Setup
sc = basicSolarSail()
date = "2023-03-21T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-22T12:30:00"
endEpoch = utc2et(endDate)  # end date
eph = twoBodyEarthEphemeride(epoch, endEpoch)
ν = get_heliocentric_position(eph, epoch)

# Parameter Setup
X0 = @MArray [42164.0, 0.65, 2*pi/180, 0.0, 0.0, 0.0] # Initial
XT = @SVector [42764.0, 0.7, 20*pi/180, 270.0*pi/180, 90*pi/180] # Targets # note that targets has 5 elements, while X0 has 6

# initialize parameters
global params = QLawParams(sc, eph, eph.t0, XT[1], XT[2], XT[3], XT[4], XT[5])
x, dx = calculate_Q_derivative(X0)

unload("naif0012.tls")
unload("de440.bsp")