# TEST SCRIPT FOR USING FRAMETRANSFORMATIONS.jl
using DrWatson
using BenchmarkTools
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
using FrameTransformations

## SPICE SETUP
furnsh(datadir("naif0012.tls"))
furnsh(datadir("de440.bsp"))
## END SPICE SETUP

# Note: SPICE treats ICRF and the J2000 frame as equivalent
epoch = 0 # start date in seconds past J2000

# Initialize some things
Orient.init_eop(datadir("iau2000a.eop.dat"))
FS = FrameSystem{4, Float64}()

# Define all axes and points needed using Spice ID's
@axes GCRF 1 GeocentricCelestialReferenceFrame
@axes ITRF 6 InternationalTerrestrialReferenceFrame
@point Earth 399
@point Spacecraft -1_900_000

# Start adding axes and points to the FrameSystem container
add_axes_inertial!(FS, GCRF)
add_point_root!(FS, Earth, GCRF)
add_axes_itrf!(FS, ITRF, GCRF) # adding the Earth-fixed reference
add_point_updatable!(FS, Spacecraft, Earth, GCRF) # adding the spacecraft as an updateable point in inertial coords

# Create our spacecraft state (Initial state from case C in q-law test scripts)
x = [-3562.1852482269273, 7107.273475624424, 0.0, -8.253230300372966, -4.136541998424384, 0.09232823522017061]
update_point!(FS, Spacecraft, x, epoch)

# Rotate this into ITRF
vector6(FS, Earth, Spacecraft, ITRF, epoch) # syntax is vector6(FrameSystem, parent, pointToExpress, frameToExpressWithin, secondsPastJ2000)