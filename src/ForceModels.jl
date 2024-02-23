#= pseudocode for earth force model
Source: Normalization of Gravitational Acceleration Models (Eckman, Brown, Adamo, 2011)

General input:
rb = [xb, yb, zb] = central-body-FIXED coordinates

Define:
r = norm(rb)

case check and trig variable assignment:
if xb=yb=0;
    sinθ = 0
    cosθ = 1
else
    sinθ = yb/sqrt(xb^2+yb^2)
    cosθ = xb/sqrt(xb^2+yb^2)
end

if r = 0
    sinΦ = 0
    cosΦ = 1
else
    sinΦ = zb/r
    cosΦ = sqrt(xb^2+yb^2) / r
end

Potential Function, V:
V = 

=#