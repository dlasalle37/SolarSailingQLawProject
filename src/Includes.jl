# Include this file in a script to include all dependencies below

# Insert all dependencies here:
# Usings
using DifferentialEquations
using Plots
using SPICE
using StaticArrays
using FiniteDiff
using ForwardDiff
using Statistics
using LinearAlgebra
using DelimitedFiles
using Infiltrator
using DiffEqCallbacks
using SparseArrays: SparseMatrixCSC, spzeros
using FrameTransformations, IERSConventions

#Includes
include("constants.jl")
include("splines.jl")
include("basicSolarSail.jl")
include("Ephemerides.jl")
include("utils.jl")
include("NormalizedGravityModel.jl")
include("QLawParams.jl")
include("EOMs.jl")
include("QLaw.jl")
include("QLawIntegrator.jl")
include("LineSearch.jl")

# Setups
