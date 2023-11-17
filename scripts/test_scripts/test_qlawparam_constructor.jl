using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP

# Spacecraft/eph Setup
sc = basicSolarSail()
date = "2023-03-21T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-22T12:30:00"
endEpoch = utc2et(endDate)  # end date
eph = twoBodyEarthEphemeride(epoch, endEpoch)

# Parameter Setup
X0 = [42164.0, 0.65, 2*pi/180, 0.0, 0.0, 0.0] # Initial
XT = [42764.0, 0.7, 20*pi/180, 270.0*pi/180, 90*pi/180] # Targets # note that targets has 5 elements, while X0 has 6
oetols = [10, 0.01, 0.1, 0.1, 0.1]

# initialize parameters
params = QLawParams(sc, eph, eph.t0, X0, XT, oetols)