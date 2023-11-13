using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

x = [2.6; 1.4; 9.8; 7.1; 0.1; 1.2] # state vector
global params = [2]
function Qf(x)

    y = sqrt(x[6])+x[5]*cos(x[4])
    
    Q = 3*x[1]+4*x[2]/(log(x[3])*y)-1
    return Q
end
FiniteDiff.finite_difference_gradient(Qf, x)