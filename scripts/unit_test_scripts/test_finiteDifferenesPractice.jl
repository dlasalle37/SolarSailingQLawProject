using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))
function Qf(x, params)

    y = sqrt(x[6])+x[5]*cos(x[4])*params[1]
    
    Q = 3*x[1]+4*x[2]/(log(x[3])*y)-1
    return Q
end

x = @MArray [2.6; 1.4; 9.8; 7.1; 0.1; 1.2] # state vector
params = [2]

R = FiniteDiff.finite_difference_gradient(x->Qf(x, params), x)