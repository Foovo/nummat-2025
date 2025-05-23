module integrali

using QuadGK #for calculating actual value of integral

export Interval, gauss_legendre, compute_integral_gl, find_best_split

struct Interval
    min
    max
end

"""
Transform original function to integrate it on an arbitrary interval
"""
function transform(f::Function, a, b)
    return t -> f(((b - a) * t + (b + a)) / 2) * (b - a) / 2
end

function gauss_legendre(int::Interval, f::Function)
    a = int.min
    b = int.max

    g = transform(f, a, b)

    return g(-sqrt(3)/3) + g(sqrt(3)/3)
end

"""
Compute the integral of a function `f` on a given interval `int` using Gauss Legendre 
by splitting the interval into `steps` steps.
"""
function compute_integral_gl(f::Function, int::Interval; steps=10)

    h = (int.max - int.min)/steps

    int_start = int.min
    gl_res = 0

    for _ = 0:steps-1
        gl_res += gauss_legendre(Interval(int_start, int_start+h), f)
        int_start += h 
    end

    return gl_res
end

function find_best_split(f::Function, int::Interval; maxiters=1000)

    actual_res, _ = quadgk(f, int.min, int.max)

    steps = 1
    while steps < maxiters
        gl_res = compute_integral_gl(f, int, steps=steps)
        steps += 1

        if isapprox(actual_res, gl_res; atol=1e-10)
            return gl_res, steps, abs(actual_res - gl_res)
        end
    end

    throw("Couldn't approximate the integral well enough")
    
end

end # module integrali
