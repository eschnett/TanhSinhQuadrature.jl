module TanhSinhQuadrature

using LinearAlgebra: norm

export Quadrature, quadts

# ∫dx f(x) = ∫dt f(t) dx/dt
#
# x(t) = tanh π/2 sinh t
# ẋ(t) = π/2 cosh t / (cosh π/2 sinh t)^2
#
# `sinh` grows exponentially.
# `tanh` limits to the range (-1, +1).

# An integration point, defined by its coordinate `x` and weight `w`
struct Point{T}
    x::T
    w::T
end

# We integrate in levels. Each successivel level refines the previous
# one, and only contains the additional points and weights.
struct Level{T}
    points::Vector{Point{T}}
end

# A complete quadrature scheme. `h` is the base step size for level 0.
struct Quadrature{T}
    h::T
    levels::Vector{Level{T}}
end

@inline ordinate(t::T) where {T<:Real} = tanh(T(π) / 2 * sinh(t))
@inline weight(t::T) where {T<:Real} = T(π) / 2 * cosh(t) / cosh(T(π) / 2 * sinh(t))^2

# Find the step size `h` at which a single step suffices to exhaust
# the numerical precision
function find_h(T::Type)
    i = 1
    while true
        t = T(i)
        w = weight(t)
        w < eps(T) && break
        i += 1
    end
    return T(max(1, i - 1))
end

function Level(h::T) where {T<:Real}
    points = Point{T}[]
    i = 1
    while true
        t = i * h
        x = ordinate(t)
        w = weight(t)
        i > 0 && abs(x) == 1 && break
        i > 0 && w == 0 && break
        push!(points, Point{T}(x, w))
        i += 2
    end
    reverse!(points)            # improve summation accuracy
    return Level{T}(points)
end

function Quadrature{T}(nlevels::Int=20) where {T<:Real}
    h = find_h(T)
    levels = Level{T}[]
    for level in 1:nlevels
        push!(levels, Level(h / 2^level))
    end
    return Quadrature{T}(h, levels)
end

function quadts(f, xmin::T, xmax::T, quad::Quadrature{T}; atol::T=zero(T),
                rtol::T=atol > 0 ? zero(T) : sqrt(eps(T))) where {T<:Real}
    h = quad.h

    x = (xmin + xmax) / 2
    w = weight(x)
    s = h * w * f(x)
    levels = 0
    error = T(Inf)

    for level in quad.levels
        h /= 2
        levels += 1
        sold = s

        s = zero(T)
        for p in level.points
            xm = xmin + (xmax - xmin) / 2 * (1 - p.x)
            xp = xmax - (xmax - xmin) / 2 * (1 - p.x)
            w = p.w
            s += w * (f(xm) + f(xp))
        end
        s = h * s + sold / 2

        tol = max(norm(s) * rtol, atol)
        error = norm(s - sold)
        levels ≥ 4 && error ≤ tol && break
    end

    return (result=s, error=error, levels=levels)
end

# function quadts(f, xmin::NTuple{D,T}, xmax::NTuple{D,T}; h::NTuple{D,T},
#                 kmax::CartesianIndex{D}) where {D,T<:Real}
#     s = zero(T)
#     for k in (-kmax):kmax
#         x = ntuple(d -> tanh(T(π) / 2 * sinh(h[d] * Tuple(k)[d])), D)
#         w = prod(ntuple(d -> T(π) / 2 * h[d] * cosh(h[d] * Tuple(k)[d]) / cosh(T(π) / 2 * sinh(h[d] * Tuple(k)[d]))^2, D))
#         s += w * f(x...)
#     end
#     return s
# end
# 
# function quadts(f, xmin::NTuple{D,T}, xmax::NTuple{D,T}; h::T, kmax::Int) where {D,T<:Real}
#     return quadts(f, xmin, xmax; h=ntuple(d -> h, D),
#                   kmax=CartesianIndex(ntuple(d -> kmax, D)))
# end

# function quadts(f, xmin::NTuple{D,T}, xmax::NTuple{D,T}; levels::Integer=4) where {D,T<:Real}
#     return quadts(f, xmin, xmax; levels=Int(levels))
# end

end