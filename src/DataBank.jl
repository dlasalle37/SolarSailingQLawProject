"""
A struct to store propagation data

    TO-DO:
        -Add this struct into QLawParams initialize in QLawParamsConstructor
        -Within QLawEOM, push each sail angle, and time to struct for saving to file later
            - will need to figure out how to only do this at the saveat times so the sail angles correspond to solution states
              if this isnt possible or hard, just remove the databank and figure out a better way to do write the angles to a file
"""
mutable struct SolarSailDataBank
    t::Vector{Vector{Float64}}
    sail_angles::Vector{Vector{Float64}}
end

# Constructor
function SolarSailDataBank()
    t = Vector{Vector{Float64}}(undef, 0)
    sail_angles = Vector{Vector{Float64}}(undef, 0)
    SolarSailDataBank(t, sail_angles)
end