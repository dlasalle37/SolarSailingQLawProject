
"""
    CubicSpline

A cubic spline interpolant.

# Fields
- `xs::AbstractVector`: The independant variable values
- `y::Union{AbstractArray, Nothing}`: An intermediate vector for storing interpolation result
- `c::AbstractArray`: The coefficients of the interpolant
"""
struct CubicSpline{
    T <: AbstractVector,
    U <: Union{AbstractArray, Nothing},
    V <: AbstractArray,
}

    # Cubic spline independant variables
    xs::T

    # Intermediate vector for storing interpolation result
    y::U

    # Coefficients
    c::V
end

"""
    CubicSpline(
        xs::AbstractVector,
        ys::Union{Vector{AbstractMatrix, VectorOfArray, AbstractArray}}
        dy0::AbstractArray,
        dyf::AbstractArray,
    )

Constructs a cubic spline interpolant.

# Arguments
- `xs::AbstractVector`: The independant variable values
- `ys::Union{Vector{AbstractMatrix, VectorOfArray, AbstractArray}}`: An array of dependant variable values, where each column is a different variable that will be interpolated.
- `dy0::AbstractVector`: The first derivative of the dependant variables at the first point
- `dyf::AbstractVector`: The first derivative of the dependant variables at the last point

# Returns
- `CubicSpline`: A cubic spline interpolant

# Throws
- `ArgumentError`: If `ys` does not have the same number of rows as the length of `xs`
- `ArgumentError`: If `ys` does not have the same number of columns as the length of `dy0`
- `ArgumentError`: If `ys` does not have the same number of columns as the length of `dyf`
"""
function CubicSpline(
    xs::AbstractVector,
    ys::Vector{T},
    dy0::AbstractVector,
    dyf::AbstractVector,
) where {T <: AbstractVector}
    return CubicSpline(xs, transpose(VectorOfArray(ys)), dy0, dyf)
end

function CubicSpline(
    xs::AbstractVector,
    ys::Union{AbstractMatrix, VectorOfArray},
    dy0::AbstractVector,
    dyf::AbstractVector,
)
    # Get number of points
    nPoinxs = length(xs)

    # Get shape of ys
    ys_shape = size(ys)

    # Check argumenxs
    if nPoinxs != ys_shape[1]
        throw(ArgumentError("ys should have the same number of rows as the length of xs"))
    end
    if length(dy0) != ys_shape[2]
        throw(ArgumentError("ys should have the same number of columns as the length of v0"))
    end
    if length(dyf) != ys_shape[2]
        throw(ArgumentError("ys should have the same number of columns as the length of vf"))
    end

    # Allocate matrix of coefficienxs
    n       = nPoinxs - 1
    c       = zeros(4*n, ys_shape[2])

    # Allocate A and b matricies
    A       = spzeros(4*n, 4*n)
    b       = zeros(4*n)

    # Loop through position positions and compute splines
    @inbounds for i in axes(ys, 2)
        # Fill A and b
        _fill_Ab!(A, b, xs, view(ys, :, i), dy0[i], dyf[i])

        # Set coefficienxs
        c[:,i] .= A \ b
    end
    return CubicSpline(xs, Vector{Float64}(undef, ys_shape[2]), c)
end

"""
    CubicSpline(
        xs::AbstractVector,
        ys::Union{Vector{AbstractVector, VectorOfArray, AbstractArray}}
        dy0::AbstractVector,
        dyf::AbstractVector,
    )

Constructs a cubic spline interpolant specialized for the case of a single dependant variable.

# Arguments
- `xs::AbstractArray`: An array of independant variable values
- `ys::Union{Vector{AbstractMatrix, VectorOfArray, AbstractArray}}`: An array of dependant variable values, where each column is a different variable that will be interpolated.
- `dy0::AbstractArray`: The first derivative of the dependant variables at the first point
- `dyf::AbstractArray`: The first derivative of the dependant variables at the last point

# Returns
- `CubicSpline`: A cubic spline interpolant

# Throws
- `ArgumentError`: If `ys` is not the same length as `xs`
"""
function CubicSpline(
    xs::AbstractVector,
    ys::AbstractVector,
    dy0::Real,
    dyf::Real,
)
    # Get number of points
    nPoinxs = length(xs)

    # Check argumenxs
    if nPoinxs != length(ys)
        throw(ArgumentError("ys and xs should be the same length."))
    end

    # Allocate matrix of coefficienxs
    n       = nPoinxs - 1
    c       = zeros(4*n)

    # Allocate A and b matricies
    A       = spzeros(4*n, 4*n)
    b       = zeros(4*n)

    # Fill A and b
    _fill_Ab!(A, b, xs, ys, dy0, dyf)

    # Set coefficienxs
    c .= A \ b

    return CubicSpline(xs, nothing, c)
end

"""
    _fill_Ab!(A, b, xs, ys, dy0, dyf)

Fills the A and b matricies in-place for the cubic spline interpolation problem.

# Arguments
- `A::AbstractMatrix`: The A matrix
- `b::AbstractVector`: The b vector
- `xs::AbstractVector`: The independant variable values
- `ys::AbstractVector`: The dependant variable values
- `dy0::Real`: The first derivative of the dependant variables at the first point
- `dyf::Real`: The first derivative of the dependant variables at the last point

# Returns
- `nothing`
"""
function _fill_Ab!(A, b, xs, ys, dy0, dyf)
    # Get n
    n  = length(xs) - 1

    # Fill A and B with 0th order constrainxs
    A .= 0.0
    b .= 0.0
    for j in 1:n
        x1 = xs[j]
        x2 = xs[j + 1]
        y1 = ys[j]
        y2 = ys[j + 1]

        A[1 + 2*(j - 1), 1 + 4*(j - 1)] = x1*x1*x1
        A[1 + 2*(j - 1), 2 + 4*(j - 1)] = x1*x1
        A[1 + 2*(j - 1), 3 + 4*(j - 1)] = x1
        A[1 + 2*(j - 1), 4 + 4*(j - 1)] = 1.0
        A[2 + 2*(j - 1), 1 + 4*(j - 1)] = x2*x2*x2
        A[2 + 2*(j - 1), 2 + 4*(j - 1)] = x2*x2
        A[2 + 2*(j - 1), 3 + 4*(j - 1)] = x2
        A[2 + 2*(j - 1), 4 + 4*(j - 1)] = 1.0
        b[1 + 2*(j - 1)]                = y1
        b[2 + 2*(j - 1)]                = y2
    end

    # Fill A and B with 1st order constrainxs
    for j in 1:n-1
        x  = xs[j + 1]
        A[2*n + j, 1 + 4*(j - 1)]           = 3.0*x*x
        A[2*n + j, 2 + 4*(j - 1)]           = 2.0*x
        A[2*n + j, 3 + 4*(j - 1)]           = 1.0
        A[2*n + j, 5 + 4*(j - 1)]           = -3.0*x*x
        A[2*n + j, 6 + 4*(j - 1)]           = -2.0*x
        A[2*n + j, 7 + 4*(j - 1)]           = -1.0
    end

    # Fill A and B with 2nd order constrianxs
    for j in 1:n-1
        x  = xs[j + 1]
        A[2*n + n - 1 + j, 1 + 4*(j - 1)]           = 6.0*x
        A[2*n + n - 1 + j, 2 + 4*(j - 1)]           = 2.0
        A[2*n + n - 1 + j, 5 + 4*(j - 1)]           = -6.0*x
        A[2*n + n - 1 + j, 6 + 4*(j - 1)]           = -2.0
    end

    # Add final constraint to force first and last 1st deriv to equal velocity
    x1  = xs[1]
    x2  = xs[end]

    A[end - 1, 1]      = 3.0*x1*x1
    A[end - 1, 2]      = 2.0*x1
    A[end - 1, 3]      = 1.0
    A[end, end - 3]    = 3.0*x2*x2
    A[end, end - 2]    = 2.0*x2
    A[end, end - 1]    = 1.0
    b[end - 1]         = dy0
    b[end]             = dyf

    return nothing
end

"""
    interpolate(spline::CubicSpline, x)

Interpolates the spline at the given value of independant variable `x`.

# Arguments
- `spline::CubicSpline`: The spline to interpolate
- `x::Real`: The value of the independant variable to interpolate to

# Returns
- `SVector`: The interpolated values of the dependant variables
"""
function interpolate(spline::CubicSpline{T,U,V}, x) where {T,U,V}
    # Check that x is in bounds
    if x < spline.xs[1] || x > spline.xs[end]
        throw(DomainError("x is outside of the spline's domain."))
    end

    # Find relevent polynomial index
    idx   = 0
    @inbounds for i in 1:length(spline.xs) - 1
        if x >= spline.xs[i] && x <= spline.xs[i + 1]
            idx = i
            break
        end
    end
    idx == 0 && throw(ErrorException("Could not find bracketing interval in xs for x."))

    # Compute interpolants
    xx   = x*x
    xxx  = xx*x
    xvec = SVector(xxx, xx, x, 1.0)
    cmat = transpose(view(spline.c, 1 + 4*(idx - 1):4*idx, :))
    mul!(spline.y, cmat, xvec)
    @infiltrate false
    return SA[spline.y[1]; spline.y[2]; spline.y[3]] 
    #return SVector(spline.y...)
end

"""
    interpolate(spline::CubicSpline, x::ForwardDiff.Dual)

Interpolates the spline at the given value of independant variable `x`, where x is a ForwardDiff Dual type.

# Arguments
- `spline::CubicSpline`: The spline to interpolate
- `x::ForwardDiff.Dual`: The value of the independant variable to interpolate to

# Returns
- `SVector`: The interpolated values of the dependant variables
"""
function interpolate(spline::CubicSpline{T,U,V}, x::ForwardDiff.Dual) where {T,U,V}
    # Check that x is in bounds
    if x < spline.xs[1] || x > spline.xs[end]
        throw(DomainError("x is outside of the spline's domain."))
    end

    # Find relevent polynomial index
    idx   = 0
    @inbounds for i in 1:length(spline.xs) - 1
        if x >= spline.xs[i] && x <= spline.xs[i + 1]
            idx = i
            break
        end
    end
    idx == 0 && throw(ErrorException("Could not find bracketing interval in xs for x."))

    # Compute interpolants
    xx   = x*x
    xxx  = xx*x
    xvec = SVector(xxx, xx, x, 1.0)
    cmat = transpose(view(spline.c, 1 + 4*(idx - 1):4*idx, :))
    y    = Vector{typeof(x)}(undef, size(cmat, 1))
    mul!(y, cmat, xvec)
    @infiltrate false
    return SVector(y...)
end


"""
    interpolate(spline::CubicSpline, x)

Interpolates the spline at the given value of independant variable `x`.

This method is specialized for the case where the spline only has one dependant variable.

# Arguments
- `spline::CubicSpline{T,Nothing,V}`: The spline to interpolate
- `x::Real`: The value of the independant variable to interpolate to

# Returns
- `Real`: The interpolated values of the dependant variables
"""
function interpolate(spline::CubicSpline{T,U,V}, x) where {T,U <: Nothing,V}
    # Check that x is in bounds
    if x < spline.xs[1] || x > spline.xs[end]
        throw(DomainError("x is outside of the spline's domain."))
    end

    # Find relevent polynomial index
    idx   = 0
    @inbounds for i in 1:length(spline.xs) - 1
        if x >= spline.xs[i] && x <= spline.xs[i + 1]
            idx = i
            break
        end
    end
    idx == 0 && throw(ErrorException("Could not find bracketing interval in xs for x."))

    # Compute interpolants
    xx   = x*x
    xxx  = xx*x
    xvec = SVector(xxx, xx, x, 1.0)
    cvec = view(spline.c, 1 + 4*(idx - 1):4*idx)
    return dot(cvec, xvec)
end

"""
    function getPositionPartials(spline::CubicSpline, x)
Retrieves the partials drdx, dydx, dzdx of the CubicSpline.

INPUTS: spline: a CubicSpline
        t: time in SPICE ephemeris time to get state from

OUTPUTS: 
        - 3-element vector where each element is a coordinate direction's partial wrt the independent variable x. 
"""
function getPositionPartials(spline::CubicSpline, x)

    # Check that t is in [0, 1]
    if x < 0.0 || x > 1.0 # Not working with Symbolics for sparsity detection atm
        throw(DomainError("Dependent variable is outside of CubicSpline bounds."))
    end

    # Find relevent polynomial index
    idx   = 0
    @inbounds for i in 1:length(spline.xs) - 1
        if x >= spline.xs[i] && x <= spline.xs[i + 1]
            idx = i
            break
        end
    end
    idx == 0 && throw(ErrorException("Could not find bracketing interval in xs for x."))  

    xx   = x*x

    # Compute Partials
    xx   = x*x
    # For a cubic ax^3+bx^2+cx+d, the derivative takes the form 3ax^2+2bx+c
    xvec = SVector(3*xx,2*x, 1.0)
    cmat = transpose(view(spline.c, 1 + 4*(idx - 1):4*idx-1, :)) #only taking (a, b, c) for each component 
    drdx = Vector{typeof(x)}(undef, size(cmat, 1))
    mul!(drdx, cmat, xvec)

    return SA[drdx[1], drdx[2], drdx[3]]
    #return SVector(drdx...)

end
