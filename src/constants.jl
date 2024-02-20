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
"""
Dictionary containing solar flux constant (G0) at certain celestial bodies (units: kgkm/s^2)
Example calculation for earth: from Oguri: Solar Sailing Primer Vector Theory
    G0 = E_e*d_e^2/c
    E_e: 1366 kg/s^2 is the solar constant at Earth (commonly given in W/m^2)
    d_e: 1.496E8 km is the sun-Earth distance
    c = 2.998E5 km/s is the speed of light
    G0_e = 1366*(1.496E8)^2/2.998E5 = 1.02E14 kgkm/s^2
"""
global SOLAR_FLUXES=Dict{Int, Any}(  # Source: NASA solar system fact sheets https://nssdc.gsfc.nasa.gov/planetary/factsheet
    10 => nothing, 
    399 => 1.02E14, #[kgkm/s^2]
    499 => 586.2*(1.52*1.496E8)^2/2.998E5, #mars: E_m=586.2W/m^2, d_m = 1.52AU, c = 2.998E5
)

global CELESTIAL_RADII=Dict{Int, Float64}( # Source: NASA solar system fact sheets https://nssdc.gsfc.nasa.gov/planetary/factsheet
    10 => 695700, #Sun [km]
    399 => 6378.0, #Earth
    301 => 1737.4, #Moon
    499 => 3389.5 #mars
)