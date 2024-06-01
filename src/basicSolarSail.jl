# Type declarations

abstract type SRPAccelerationModel        end # Supertype for acceleration models
abstract type FPOF <:SRPAccelerationModel end # Subtype for flat-plate optical force model

# The basicSolarSail struct
struct basicSolarSail
    areaParam::Float64          # A/m [km^2/kg]
    C::Vector{Float64}          # contains all 3 (C1, C2, C3) condensed parameters (Default is NASA NEA Scout)
    accelerationModel::Symbol   # acceleration model to be used; default is FPOF for flat-plate optical force (Model used in Oguri)
end

# constructor
function basicSolarSail(; areaParam=100*5.4E-6, C1=1.711, C2=0.002, C3=0.145, accelerationModel=:FPOF)
    if (C1==0) && (C2==0)
        error("Check sail parameters, C1 and C2 cannot both be zero")
    elseif !(0<=C1<=2) || !(0<=C2<=1)
        error("Values out of bounds. C1∈[0,2] and C2∈[0,1]")
    end

    condensedParams = [C1,C2,C3]

    # Create struct
    basicSolarSail(areaParam, condensedParams, accelerationModel)
end