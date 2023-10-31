using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

date = "2023-03-01T12:30:00" 
epoch = utc2et(date)  # start date
endDate = "2023-03-01T12:30:00"
endEpoch = utc2et(endDate)  # end date
a = twoBodyEarthEphemeride(epoch, endEpoch)

ν = earth_heliocentric_position(a, endEpoch)