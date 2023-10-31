# Include this file in a script to include all dependencies below

# Insert all dependencies here:
# Usings
using DifferentialEquations
using Plots
using SPICE

#Includes
include("constants.jl")
include("Ephemerides.jl")
include("basicSolarSail.jl")
include("utils.jl")



# Setups
## SPICE SETUP
furnsh("naif0012.tls")
furnsh("de440.bsp")
## END SPICE SETUP