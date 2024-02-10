# Planetary constants
global EARTH_MU = 398600.4418 # km^3/s^2
global EARTH_SOLAR_FLUX = 1.02E14 
# Approximate constants for Earth's heliocentric orbit
global SUN_MU = 1.327E11  #km^3/s^2
global EARTH_HELIOC_ECC = 0.0167
global EARTH_HELIOC_SEMIMAJOR = 149598023.0  # km
global SUN_EARTH_DIST_NOMINAL = 149597870.691 # km

"""
Quick dictionary for the gravitational parameters of various celestial bodies, indexed by their SPICE id's (all in km^3/s^2)
    - Note: Maybe make this dictonary more sophisticated by having each entry be a struct containing a bunch of planetary data rather
            than just the gravitational parameter
    - from: https://nssdc.gsfc.nasa.gov/planetary/factsheet
"""
global GRAV_PARAMS=Dict{Int, Float64}(
    10 => 1.32712440018E11, #Sun
    399 => 398600.4418, #Earth
    301 => 4902.800118, #moon
    499 => 4.28284E4 #Mars
)

global SOLAR_FLUXES=Dict{Int, Any}(
    10 => nothing, 
    399 => 1.02E14 #[kgkm/s^2]
)