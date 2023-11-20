struct basicSolarSail
    areaParam::Float64
    C::Vector{Float64} # contains all 3 (C1, C2, C3) condensed parameters (Default is NASA Nea Scout)
    accelerationModel::Symbol  # acceleration model to be used; default is FPOF for flat-plate optical force (Model used in Oguri)
    method::Symbol  # which states will be used? (Oguri for lambda, Kep for regular RAAN)
end

# constructor
function basicSolarSail(; areaParam=100*5.4E-6, C1=1.711, C2=0.002, C3=0.145, accelerationModel=:FPOF, method=:Oguri)
    if (C1==0) && (C2==0)
        error("Check sail parameters, C1 and C2 cannot both be zero")
    elseif !(0<=C1<=2) || !(0<=C2<=1)
        error("Values out of bounds. C1∈[0,2] and C2∈[0,1]")
    end

    condensedParams = [C1,C2,C3]

    # Create struct
    basicSolarSail(areaParam, condensedParams, accelerationModel, method)
end

"""
Depreciated
"""
function calculate_sail_constants(ρ, s, ϵf, ϵb, Bf, Bb)
    #= 
    This function executes in the constructor.
    Uses the sail parameters to calculate the condensed parameters C1, C2, C3
    Formulation discussed in Solar Sailing Primer Vector Theory: Indirect Trajectory
    Optimization with Practical Mission Considerations by Oguri, McMahon, Lantoine.
    =#

    C1 = 2*ρ*s
    C3 = 1-ρ*s
    C2 = Bf*ρ*(1-s)+(1-ρ)*(Bf*ϵf-Bb*ϵb)/(ϵf+ϵb)

    if (C1==0) && (C2==0)
        error("Check sail parameters, C1 and C2 cannot both be zero")
    elseif !(0<=C1<=2) || !(0<=C2<=1)
        error("Values out of bounds. C1∈[0,2] and C2∈[0,1]")
    end
    return C1, C2, C3
end


