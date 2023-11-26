using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

function Qf(x, params) # Dummy function

    y = sqrt(x[6])+x[5]*cos(x[4])*params[1]
    
    Q = 3*x[1]+4*x[2]/(log(x[3])*y)-1
    return Q
end
params = [2, 3]
x = @MArray [2.6; 1.4; 9.8; 7.1; 0.1; 1.2] # state vector

# ForwardDiff specific stuff
out = similar(x)
cfg = ForwardDiff.GradientConfig(Qf, x) # GradientConfig

# Solving
soln = ForwardDiff.gradient!(out, x->Qf(x, params), x, cfg)