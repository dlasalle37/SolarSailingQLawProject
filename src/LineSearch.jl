"""
    gss(f, a, b; tol, maxiter)
Find local minimum of function f using the Golden Section Search method

# Inputs
    f: univariate function handle
    a: lower search bound
    b: upper search bound
    tol: specified tolerance
    maxiter: maximum iterations before quit

# Returns
    xstar: minimizing value of function f
    f(xstar): local minimum value of function f

# Source
    https://indrag49.github.io/Numerical-Optimization/solving-one-dimensional-optimization-problems.html#golden-section-search-method
"""
function gss(f::Function, a, b; tol=1e-5, maxiter=500)
    gr = (sqrt(5) - 1)/2
    i = 0
    x1 = gr*a + (1-gr)*b
    x2 = (1-gr)*a + gr*b
    fx1 = f(x1)
    fx2 = f(x2)
    while b - a > tol
        
        if fx1 < fx2
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = gr*a + (1-gr)*b
            fx1 = f(x1)
        else
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = (1-gr)*a + gr*b
            fx2 = f(x2)
        end
        if i >= maxiter # break if not converged
            return NaN, NaN, i
            break
        end
        i += 1

    end

    xstar = (a+b)/2
    fxstar = f(xstar)
    return xstar, fxstar, i
end