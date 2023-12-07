using DrWatson
@quickactivate "SolarSailingQLawProject"
include(srcdir("Includes.jl"))

sc = basicSolarSail() # [C1; C2; C3] = [1.711, 0.002, 0.145]
p = [1.0; 1.5; 2.0] # arbitrary P
alphastar = calculate_alpha_star(p, sc)
println("Function result: $alphastar")

# Check calc
C1 = sc.C[1]; C2 = sc.C[2]; C3 = sc.C[3]
k = p[1]/sqrt(p[2]^2+p[3]^2)
α0 = atan(1/4 * (-3k + sqrt(8+9*k^2)))
a = α0
F_a = C1*cos(a)*(3*sin(a)^2 - 1) - C2*cos(2*a) + k*sin(a)*(C3 + 2*C2*cos(a) + 3*C1*cos(a)^2)
F_a_a = 9*C1*k*cos(a)^3 - 2*C1*sin(a) - 2*C2*k + 4*C2*k*cos(a)^2 + 4*C2*cos(a)*sin(a) + 9*C1*cos(a)^2*sin(a) - 6*C1*k*cos(a) + C3*k*cos(a)
res = α0 - F_a/F_a_a
println("Recalculated value from symbolic derivatives: $res")
assessment = res==alphastar
println("Test result (true for pass): $assessment")