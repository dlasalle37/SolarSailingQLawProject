include("Ephemerides.jl")


"""
    third_body_perturb(eph::Ephemeride, perturbingBody, et)
Calculate the perturbing acceleration on a spacecraft due to a third body (e.g. Sun, Moon, etc.)

# Inputs:
    pos: spacecraft position vector around in an inertial frame centered at the body it is orbiting (e.g. J2000 for Earth-centered sc)
    et: ephemeris time to get acceleration vector at
    eph::Ephemeride: Struct carrying the information about the third body (eph.frame must match frame of pos)
    perturbingBody: NAIF/SPICE ID of third body (e.g. 10 for sun, 399 for Earth, etc.)

# Outputs:
    a_tb: perturbing acceleration due to third body (km/s^2)
"""
function third_body_perturb(pos, et, eph::Ephemeride, perturbingBody)
    negate = false # bool for determining if it is necessary to negate the position vector returned from the ephemeride
    # Make sure that the ephemeride contains info about the desired perturbing body
    if perturbingBody != eph.targ   #Usually, the perturbing body would be  the ephemeride target, but note if it is not
        if perturbingBody == eph.obs
            # If the observer is the desired third body, then the relative position vectors will need to be 
            # negated for the vector math to work
            negate = true           
        else
            throw(ArgumentError("State information for celestial body $perturbingBody not found in Ephemeride struct"))
        end
    end

    # Get relative position vectors
    r_sc_cb = pos; # making this clear, r_sc_cb is position of spacecraft rel. to central body (assumed inertial in cb-centered frame)

    r_tb_cb = interpolate(eph, et); # position vector pointing from eph.obs to eph.targ
    if negate; r_tb_cb = -1*r_tb_cb end # negate if necessary so it points from central body to third body

    r_cb_tb = -r_tb_cb

    r_sc_tb = r_sc_cb + r_cb_tb;
    r = norm(r_sc_tb)

    # calculate accelerative force
    mu_tb = get_gm(perturbingBody) # gravitational parameter of third body
    a_tb = mu_tb/r^3 * r_sc_tb 

    return a_tb


end