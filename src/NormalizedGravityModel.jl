#= Gravitational Force Modeling via harmonics
Source: Fast Gravity with partials nasa technical doc, 1988, appendix A3 for code

Potential Function, U:
U = µ/r * sum[n=2 to inf](sum[m=0 to n](µ/r * (α/r)^n * P_{n,m}(x_3/r) * (C_{n,m}cos(mλ)+S_{n,m}sin(mλ)) ))

=#
struct NormalizedGravityModelData

    n::Int # max degree
    m::Int # max order

    R::Float64 # Planet average equatorial radius
    mu::Float64 # planet gravitational parameter

    C::AbstractArray # Normalized matrix of cosine coefficients
    S::AbstractArray # Normalized matrix of sine coefficients

end

"""
    NormalizedGravityModelData(order_and_degree, Cfile, Sfile, R, mu)
Constructor for Normalized Gravity model data. For now, only square models are considered (i.e same degree and order). 
This model uses normalized coefficients directly.

# Inputs:
    order_and_degree: integer value of degree and order to calculate to.
    Cfile: location of cosine coefficients
    Sfile: location of sine coefficients
    R: planetary radius [km]
    mu: planetary gravitational parameter [km^2/s^3]

# Outputs:
    NormalizedGravityModelData: struct to access all necessary data for gravity potential based calculations
"""
function NormalizedGravityModelData(order_and_degree, Cfile, Sfile; R=6738.0, mu=398600.4418)
    # Assign degree/order (assuming square for now)
    n = order_and_degree
    m = order_and_degree

    # n, m error check
    if m > n
        throw(error("Order should be less than or equal to degree"))
    end
    
    # read in coefficients
    Cfull = readdlm(Cfile, '\t', '\n'; header=false)
    Sfull = readdlm(Sfile, '\t', '\n'; header=false) 
    C = view(Cfull, 1:n, 1:m)
    S = view(Sfull, 1:n, 1:m) 

    # Create model
    NormalizedGravityModelData(
        n,
        m,
        R,
        mu,
        C,
        S
    )
end

"""
    getPot(model::NormalizedGravityModelData, x::AbstractVector, wantCentralForce::bool)
returns the gravitational acceleration vector g in earth-fixed coordinates by calculating the first partial of the gravitational potential.
"""
function getFirstPartial(
    model::NormalizedGravityModelData, # 
    x::AbstractVector, # fixed position vector
    wantCentralForce::Bool # do you want entire potential or just pertubation terms
)
    # Point to data in NormalizedGravityModelData
    Re = model.R    # planet radius
    C = model.C     # cosine coeffs.
    S = model.S     # sine coeffs
    nmax = model.n  # desired degree
    mmax = model.m  # desired order
    mu = model.mu   # planet grav. parameter


    # Calculate some initial terms
    r = norm(x)
    Xovr = x[1] / r             # x/r
    Yovr = x[2] / r             # y/r
    Zovr = x[3] / r             # z/r
    ep = Zovr                   # epsilon
    mu_over_r = mu / r          # mu/r
    mu_over_r2 = mu_over_r / r  # mu/r^2
    Reovr = Re / r              # Re/r
    Reovrn = Reovr              # (Re/r)^n (will be calculated in recursion)

    # C_tilda and S_tilda terms, equation 3-18
    Ctil = Vector{}(undef, mmax) 
    Stil = Vector{}(undef, mmax)
    Ctil[1] = 1.0; Ctil[2] = Xovr;
    Stil[1] = 0.0; Stil[2] = Yovr;

    # Initialize Summations
    if wantCentralForce; sum_init = 0; else sum_init = 1; end
    sumh = 0.0
    sumj = 0.0
    sumk = 0.0
    sumgam = sum_init

    # Generating associated legendre functions
    P = zeros(nmax, mmax) # initialize array of all legendre functions
    P[1,1] = 1.0 # seed recursion with initial values
    P[1,2] = 0.0
    P[2,1] = sqrt(3)*ep
    P[2,2] = sqrt(3)

    for n = 2:nmax-1
        nidx = n + 1; # Shift bc julia indexing starts at 1. Note anywhere where n appears in calculation, real n is used, while nidx is only for indexing
        Reovrn = Reovrn*Reovr

        # First col. terms
        P[nidx, 1] = xi(n,0)*ep*P[nidx-1, 1] - eta(n,0)*P[nidx-2, 1]

        # inner-diagonal terms
        deltan = sqrt(2*n+1) * P[nidx-1, nidx-1]    # calculating the delta here. The tech document says delta only depends on n and not the state, 
                                                    # but it depends on previous P values, which come from the state, so not sure how they are pre-allocating it
        P[nidx, nidx-1] = deltan * ep

        # Diagonal term
        P[nidx, nidx] = sqrt((2*n+1)/(2*n)) * P[nidx-1, nidx-1]

        # Initialize summations over m=1->n (m-based summations)
        sumH_n = zeta(n, 0) * P[nidx, 2] * C[nidx, 1]
        sumgam_n = P[nidx, 1] * C[nidx, 1] * (n+1)
        sumJ_n = 0.0
        sumK_n = 0.0

        if mmax > 0 # if mmax is zero, we are only using zonal coefficients, (C_1,0, C_2,0, ...) and (S_1,0, S_2,0, ....)
            for m = 1:n-2
                midx = m+1 # same deal here with index vs. actual
                # Calculate all other legendre's in row nidx
                P[nidx, midx] = xi(n, m)*ep*P[nidx-1, midx] - eta(n, m)*P[nidx-2, midx]
            end

            # Calculate new Ctilda, Stilda via recursion 
            Ctil[nidx] = Ctil[2] * Ctil[nidx-1] - Stil[2] * Stil[nidx-1]
            Stil[nidx] = Stil[2] * Ctil[nidx-1] + Ctil[2] * Stil[nidx-1]
            
            # This statement comes into play if order is less than degree. Only happens for non-square models, (e.g, (40,35), (200,100), etc.)
            if n < mmax
                lim = n
            else
                lim = mmax
            end

            # Calculate m-based summations (Equations 4-11)
            for m = 1:lim-1
                midx = m+1

                # Calculate some shorthand terms
                Pnm = P[nidx, midx]
                Cnm = C[nidx, midx]
                Snm = S[nidx, midx]
                Bnmtil = Cnm*Ctil[midx] + Snm*Stil[midx]

                # Summations (Eq 4-11 for all except sumH_n)
                sumJ_n = sumJ_n + m*Pnm*(Cnm*Ctil[midx-1] + Snm*Stil[midx-1])
                sumgam_n = sumgam_n + (n+m+1)*Pnm*Bnmtil
                sumK_n = sumK_n + m*Pnm*(Snm*Ctil[midx-1] - Cnm*Stil[midx-1])

                # Sum for Hn is given as below, note only difference from Eq4-11 is multiplication by zeta
                sumH_n = sumH_n + P[nidx, midx+1]*Bnmtil*zeta(n, m)

            end

            # add all m-based summations to n-based summations (n-based meaning summed from n=2:nmax)
            sumj = sumj + Reovrn*sumJ_n
            sumk = sumk + Reovrn*sumK_n
        end

        # These two will have values even if mmax = 0 because they are initialized at a nonzero value
        # so sum them outside of the if statement
        sumh = sumh + Reovrn*sumH_n
        sumgam = sumgam + Reovrn*sumgam_n
    end

    # Putting all the summations together into the gradient 
    Lambda = sumgam + ep*sumh
    g1 = -mu_over_r2 * (Lambda*Xovr - sumj)
    g2 = -mu_over_r2 * (Lambda*Yovr - sumk)
    g3 = -mu_over_r2 * (Lambda*Zovr - sumh)
    
    dUdx = @SVector [g1; g2; g3]

    return dUdx, P # note remove P as output later, just for testing

end

# SUPPLEMENTAL FUNCTIONS
"""
    N(n,m)
Normalization coefficient calculator for a given degree and order n,m
Note: This will easily overflow due to factorial(n+m) in denominator
# Inputs:
    - n: degree
    - m: order
# Outputs
    N: normalization factor
"""
function N(n,m)
    if m == 0
        factor = (2*n+1)^1/2
    else
        factor = sqrt( factorial(n-m)*(2n+1)*2 / factorial(n+m) )
    end
    return factor
end

"""
    normalized_legendre_generator(n, m, ε)
-recursively generate all legendre polynomials needed (each is evaluated at ε).
-Note this is currently for testing only, as legendre functions are generated within the getPot methods. This function is an isolated way to generate all
 legendre polynomials used throughout the potential partial calculations. Make sure that the matrix P in that function matches the P in this function
 for a given ε.
"""
function normalized_legendre_generator(nmax, mmax, ε)
    # init
    P = zeros(nmax, mmax)

    # Begin recursion
    P[1,1] = 1.0 #P @ n=0, m=0 (index shift for Julia)
    P[2,1] = sqrt(3)*ε
    P[2,2] = sqrt(3)
    P[1,2] = 0.0

    for n=2:nmax-1
        nidx = n + 1; # julia indexing starts at 1
        # Compute diagonal term
        P[nidx,nidx] = sqrt((2*n+1)/(2*n)) * P[nidx-1, nidx-1]

        # Compute leftmost column term
        P[nidx,1] = ε*xi(n,0)*P[nidx-1, 1] - eta(n,0)*P[nidx-2, 1]

        # Compute below-diagonal
        δ = sqrt(2*n+1) * P[nidx-1, nidx-1]
        P[nidx, nidx-1] = ε*δ

        for m = 1:n-2 #Starting at column second from the left
            midx = m+1
            ξ = xi(n,m)
            η = eta(n,m)
            P[nidx, midx] = ξ*ε*P[nidx-1, midx] - η*P[nidx-2,midx]
        end
    end


    return P
end

"""
    eta(n, m)
"""
function eta(n,m)
    η = sqrt( (n+m-1)*(2*n+1)*(n-m-1) / ((n+m)*(n-m)*(2*n-3)) )
    return η
end

"""
xi(n, m)
"""
function xi(n,m)
    ξ = sqrt( (2*n-1)*(2*n+1)/((n-m)*(n+m)) )
    return ξ
end

"""
    zeta(n,m)
parameter zeta
"""
function zeta(n, m)
    if m == 0
        ζ = sqrt( n*(n+1)/2 )
    else
        ζ = sqrt( (n-m)*(n+m+1) )
    end

    return ζ
end

"""
    upsilon(n,m)
parameter upsilon (used for second partial)
"""
function upsilon(n,m)
    if m==0
        Y = sqrt( n*(n-1)*(n+1)*(n+2)/2 )
    else
        Y = sqrt((n-m)*(n-m-1)*(n+m+1)*(n+m+2))
    end

    return Y
end