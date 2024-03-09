#= 
Gravitational Force Modeling via harmonics
Source: Fast Gravity with partials nasa technical doc, 1993. (NASA CR-188243)
=#

struct NormalizedGravityModelData

    n::Int # max degree
    m::Int # max order

    R::Float64 # Planet average equatorial radius
    mu::Float64 # planet gravitational parameter

    C::AbstractArray # Normalized matrix of cosine coefficients
    S::AbstractArray # Normalized matrix of sine coefficients

    GravityModel::Symbol # Model name (example; :EGM96)

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
    mod: Symbol for gravity model desired (defaults to EGM96)

# Outputs:
    NormalizedGravityModelData: struct to access all necessary data for gravity potential based calculations
"""
function NormalizedGravityModelData(degree, order, fileloc; R=6738.0, mu=398600.4418, GravModel=:EGM96)
    # Assign degree/order
    n = degree
    m = order
    nidx = n+1; # Accounting for julia indexing shift
    midx = m+1; # Accounting for julia indexing shift

    # read in coefficients based on model given
    Cfull, Sfull = EGM_coeff_reader(fileloc)

    # n, m error check
    if m > n
        throw(error("Order should be less than or equal to degree"))
    elseif nidx > size(Cfull, 1) || midx > size(Cfull, 2)
        throw(error("Specified degree or order larger than available coefficients"))
    end
    
    C = view(Cfull, 1:nidx, 1:midx)
    S = view(Sfull, 1:nidx, 1:midx) 

    # Create model
    NormalizedGravityModelData(
        n,
        m,
        R,
        mu,
        C,
        S,
        GravModel
    )
end

"""
    EGM_coeff_reader
File reader designed to read in and organize Cnm and Snm coefficient data for the EGM earth model. Current implementation is designed to read the specific file
    "EGM96_to360.ascii" available at "https://cddis.nasa.gov/926/egm96/getit.html"
# Inputs: 
    fileloc: string containing full local file location
# Outputs:
    C: matrix of C coefficients
    S: matrix of S coefficients
    """
function EGM_coeff_reader(fileloc)

    # Read File
    data = readdlm(fileloc)

    # Add 1 to the index columns (col. 1 and 2 for n and m, respectively) to account for julia index shift
    data[:,1:2] .+= 1

    # Instantiate C and S matrices
    nmax = Int64(data[end, 1])
    mmax = Int64(data[end, 2])
    C = Matrix{Float64}(undef, nmax, mmax)
    C[1:2,1:2] .= 0;                         #assumed zero
    S = Matrix{Float64}(undef, nmax, mmax)   # S is instantiated the same as C
    S[1:2,1:2] .= 0;  

    # In the assumed file, the C coeffs are column 3, and the stddev of the C coeffs are column 5
    # Similarly, the S coeffs are column 4 and the stddev of the S coefficients are column 6
    for row in eachrow(data)
        nidx = Int64(row[1])           # first col is n
        midx = Int64(row[2])           # second col is m
        C[nidx, midx] = row[3]  # Assign Cnm from column 3
        S[nidx, midx] = row[4]  # Assign Snm from column 4
    end
    return C, S
end

"""
    getFirstPartial(model::NormalizedGravityModelData, x::AbstractVector, wantCentralForce::bool)
returns the gravitational acceleration vector g in earth-fixed coordinates by calculating the first partial of the gravitational potential.
"""
function getFirstPartial(
    mdl::NormalizedGravityModelData, # data struct
    x::AbstractVector, # fixed position vector
    wantCentralForce::Bool # do you want entire potential or just pertubation terms
)
    # Point to data in NormalizedGravityModelData
    Re = mdl.R    # planet radius
    C = mdl.C     # cosine coeffs.
    S = mdl.S     # sine coeffs
    nmax = mdl.n  # desired degree
    mmax = mdl.m  # desired order
    mu = mdl.mu   # planet grav. parameter


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
    Ctil = Vector{}(undef, nmax) 
    Stil = Vector{}(undef, nmax)
    Ctil[1] = 1.0; Ctil[2] = Xovr;
    Stil[1] = 0.0; Stil[2] = Yovr;

    # Initialize Summations
    # gamma is the only summation that conditionally starts at 1 or 0, based on if the user wants to include the central force or not, respectively.
    if wantCentralForce; sumgam = 1.0; else sumgam = 0.0; end 
    sumh = 0.0
    sumj = 0.0
    sumk = 0.0

    # Generating associated legendre functions
    P = zeros(nmax, nmax)   # initialize array of all legendre functions. n>=m always true, so initialize as nxn to avoid errors with non-square models.
    P[1,1] = 1.0            # seed recursion with initial values
    P[1,2] = 0.0
    P[2,1] = sqrt(3)*ep
    P[2,2] = sqrt(3)

    for n = 2:nmax-1
        nidx = n + 1;           # Shift bc julia indexing starts at 1. Note anywhere where n appears in calculation, real n is used, while nidx is only for indexing
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
        # Done here because some will have a value even if mmax=0
        sumH_n = zeta(n, 0) * P[nidx, 2] * C[nidx, 1]
        sumgam_n = P[nidx, 1] * C[nidx, 1] * (n+1)
        sumJ_n = 0.0
        sumK_n = 0.0

        if mmax > 0 # if mmax is zero, we are only using zonal coefficients, (C_1,0, C_2,0, ...) and (S_1,0, S_2,0, ....)
            for m = 1:n-2
                midx = m+1 # same deal here with index vs. actual
                # Calculate all other legendre's in row nidx
                # This must be done before calculating other summations due to m+1 terms that appear in those.
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

"""
function getSecondPartial(model::NormalizedGravityModelData, x::AbstractVector, wantCentralForce::Bool)
    NOT YET WORKING
"""
function getSecondPartial(mdl::NormalizedGravityModelData, x::AbstractVector, wantCentralForce::Bool)
    # Pull data from struct
    nmax = mdl.n
    mmax = mdl.m
    Re = mdl.R
    mu = mdl.mu
    C = mdl.C
    S = mdl.S

    # Some initial/shorthand calcs
    r = norm(x)
    ep = x[3] / r
    Xovr = x[1]/r
    Yovr = x[2]/r
    Reovr = Re/r
    Reovrn = Reovr
    mu_over_r = mu/r
    mu_over_r2 = mu_over_r / r
    mu_over_r3 = mu_over_r2 / r
    # C_tilda and S_tilda terms, equation 3-18
    Ctil = Vector{}(undef, nmax) 
    Stil = Vector{}(undef, nmax)
    Ctil[1] = 1.0; Ctil[2] = Xovr;
    Stil[1] = 0.0; Stil[2] = Yovr;

    # initialize conditional n-based summations
    if wantCentralForce; sum_init=1.0; else sum_init=0.0; end
    sumgam = sum_init
    sumL = 2.0*sum_init

    # Initialize other n-based summations
    sumN = 0.0
    sumM = 0.0
    sumOmega = 0.0
    sumP = 0.0
    sumQ = 0.0
    sumR = 0.0
    sumS = 0.0
    sumT = 0.0
    sumH = 0.0

    # Initialize legendre functions organized in matrix P
    P = zeros(nmax, nmax)   # initialize array of all legendre functions. n>=m always true, so initialize as nxn to avoid errors with non-square models.
    P[1,1] = 1.0            # seed recursion with initial values
    P[1,2] = 0.0
    P[2,1] = sqrt(3)*ep
    P[2,2] = sqrt(3)
    

    for n = 2:nmax-1
        nidx = n+1              # Indexing in Julia starts at 1. Note that n is used in formulas, while nidx is for indexing only.
        Reovrn = Reovrn*Reovr   # Update (Re/r)^n
        
        # Computing first set of Legendre functions (that don't depend on m)
        # First col. terms
        P[nidx, 1] = xi(n,0)*ep*P[nidx-1, 1] - eta(n,0)*P[nidx-2, 1]

        # inner-diagonal terms
        deltan = sqrt(2*n+1) * P[nidx-1, nidx-1]    # calculating the delta here. The tech document says delta only depends on n and not the state, 
                                                    # but it depends on previous P values, which come from the state, so not sure how they are pre-allocating it
        P[nidx, nidx-1] = deltan * ep

        # Diagonal term
        P[nidx, nidx] = sqrt((2*n+1)/(2*n)) * P[nidx-1, nidx-1]

        # Initialize m-based summations
        # Sums that have a value at m=0:
        sumgam_n = P[nidx, 1] * C[nidx, 1] * (n+1)
        sumL_n = P[nidx, 1] * C[nidx, 1] * (n+1) * (n+2)
        sumP_n = P[nidx, 2] * C[nidx, 1] * zeta(n,0)
        sumM_n = P[nidx, 3] * C[nidx, 1] * upsilon(n,0)
        sumH_n = zeta(n, 0) * P[nidx, 2] * C[nidx, 1]

        # Sums that are zero at m=0:
        sumOmega_n = 0.0 
        sumN_n = 0.0
        sumQ_n = 0.0
        sumR_n = 0.0
        sumS_n = 0.0
        sumT_n = 0.0


        if mmax > 0 # skip all m-based items if we are only using zonals
            # Get rest of Legendre functions in row nidx
            for m = 1:n-2
                midx = m + 1
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
            for m = 1:lim-1
                midx = m + 1

                # Define some shorthand terms
                Pnm = P[nidx, midx]
                Cnm = C[nidx, midx]
                Snm = S[nidx, midx]
                Bnmtil = Cnm*Ctil[midx] + Snm*Stil[midx]
                
                # Summations (Eq 5-10)
                sumgam_n = sumgam_n + (n+m+1)*Pnm*Bnmtil
                sumH_n = sumH_n + P[nidx, midx+1]*Bnmtil*zeta(n, m)
                sumL_n = sumL_n + (n+m+1)*(n+m+2)*Pnm*Bnmtil
                sumP_n = sumP_n + P[nidx, midx+1]*(m+n+1)*Bnmtil * zeta(n,m)
                sumQ_n = sumQ_n + P[nidx, midx+1]*m*(Cnm*Ctil[midx-1] + Snm*Stil[midx-1]) * zeta(n,m)
                sumR_n = sumR_n + P[nidx, midx+1]*m*(Snm*Ctil[midx-1] - Cnm*Stil[midx-1]) * zeta(n,m)
                sumS_n = sumS_n + Pnm*m*(n+m+1)*(Cnm*Stil[midx-1] - Snm*Ctil[midx-1])
                sumT_n = sumT_n + Pnm*m*(n+m+1)*(Snm*Ctil[midx-1] - Cnm*Stil[midx-1])
                
                # Legendre check for sumM, because of m+2 index (avoiding bounds error) 
                if midx+2 > size(P,2)
                    sumM_n = sumM_n + 0 * Bnmtil * upsilon(n,m) # P[n,m>n] = 0 by definition, even if we haven't actually calculated out that far
                else
                    sumM_n = sumM_n + P[nidx, midx+2] * Bnmtil * upsilon(n,m)
                end

                # Ctil & Stil check
                if m >= 2 # can only calculate these terms if m>=2 (assumed Ctil[m<0] = Stil[m<0] = 0)
                    sumOmega_n = sumOmega_n + Pnm*m*(m-1)*(Cnm*Stil[midx-2] - Snm*Ctil[midx-2])
                    sumN_n = sumN_n + Pnm*m*(m-1)*(Cnm*Ctil[midx-2]+Snm*Stil[midx-2])
                end
            end

            # Add to n-based summations that were initialized at zero
            sumOmega = sumOmega + Reovrn * sumOmega_n
            sumN = sumN + Reovrn * sumN_n
            sumQ = sumQ + Reovrn * sumQ_n
            sumR = sumR + Reovrn * sumR_n
            sumS = sumS + Reovrn * sumS_n
            sumT = sumT + Reovrn * sumT_n
        end
        # Add to n-based summations that were initialized at nonzero value
        sumgam = sumgam + Reovrn * sumgam_n
        sumL = sumL + Reovrn * sumL_n
        sumP = sumP + Reovrn * sumP_n
        sumM = sumM + Reovrn * sumM_n
        sumH = sumH + Reovrn * sumH_n
    end

    # Use completed sums to obtain second partial matrix
    # Eq 5-11
    avec = @SVector [sumQ; sumR; 0]
    yvec = @SVector [sumS; sumT; 0]

    # Define some needed quantities (Recall from first partial method)
    xhat = x/r
    Lambda = sumgam + ep*sumH

    # Eq 5-14
    F = sumL + ep*(sumM*ep + 2*(sumP+sumH)) + Lambda
    G = -(sumM*ep + sumP + sumH)
    dvec = ep*avec + yvec

    # Eq 5-15
    # First added term
    FGM = @SMatrix [F G; G sumM]
    term1 = hcat(xhat, avec)*FGM*vcat(transpose(xhat), transpose(avec))

    # Second added term
    negI = @SMatrix [0 -1; -1 0] # negative transpose of I2x2
    term2 = hcat(xhat, dvec)*negI*vcat(transpose(xhat), transpose(dvec))

    # Third added term
    term3 = @SMatrix [sumN-Lambda -sumOmega sumQ; -sumOmega -(sumN+Lambda) sumR; sumQ sumR -Lambda]

    d2Udx2 = mu_over_r3 * (term1 + term2 + term3)
    return d2Udx2
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
-Note this is currently for testing only, as legendre functions are generated within the getPartial methods. This function is an isolated way to generate all
 legendre polynomials used throughout the potential partial calculations. Make sure that the matrix P in that function matches the P in this function
 for a given ε.
"""
function normalized_legendre_generator(nmax, mmax, ε)
    # init
    P = zeros(nmax, nmax)

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
Used in legendre function recursive generation
"""
function eta(n,m)
    η = sqrt( (n+m-1)*(2*n+1)*(n-m-1) / ((n+m)*(n-m)*(2*n-3)) )
    return η
end

"""
    xi(n, m)
Used in legendre function recursive generation
"""
function xi(n,m)
    ξ = sqrt( (2*n-1)*(2*n+1)/((n-m)*(n+m)) )
    return ξ
end

"""
    zeta(n,m)
parameter zeta. Used in first/second partial to calculate Hn based on normalized legendre functions and coeffs.
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