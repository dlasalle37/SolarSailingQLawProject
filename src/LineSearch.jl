"""
    gss(f, xa, xb; tol, maxiter)
Find local minimum of function f using the Golden Section Search method

# Inputs
    f: univariate function handle
    xa: lower search bound
    xb: upper search bound
    tol: specified tolerance
    maxiter: maximum iterations before quit

# Returns
    xstar: minimizing value of function f
    f(xstar): local minimum value of function f

# Source
    https://indrag49.github.io/Numerical-Optimization/solving-one-dimensional-optimization-problems.html#golden-section-search-method
"""
function gss(f, xa, xb; tol=1e-5, maxiter=500)
    gr = (sqrt(5) - 1)/2
    i = 0
    L = xb - xa
    x1 = xa + gr^2*L
    x2 = xa + gr*L
    while L > tol
        
        if f(x1) > f(x2)
            xa = x1
            x1 = x2
            L = xb - xa
            x2 = xa + gr * L
        else
            xb = x2
            x2 = x1
            L = xb - xa
            x1 = xa + gr^2 * L
        end
        if i >= maxiter # break if not converged
            return NaN, NaN
            break
        end
        i += 1

    end

    xstar = (xa+xb)/2
    fxstar = f(xstar)
    return xstar, fxstar
end