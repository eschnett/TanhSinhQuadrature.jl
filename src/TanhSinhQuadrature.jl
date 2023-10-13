module TanhSinhQuadrature

using LinearAlgebra: norm
using StaticArrays
using SimpleNonlinearSolve
################################################################################

export TSQuadrature, quadts

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
struct TSQuadrature{T}
    h::T
    levels::Vector{Level{T}}
end

@inline ordinate(t::T, p::Int=1) where {T<:Real} = tanh(T(π) / 2 * sinh(t^p))
@inline weight(t::T, p::Int=1) where {T<:Real} = (p / 2) * T(π) * t^(p - 1) * cosh(t^p) / cosh(T(π) / 2 * sinh(t^p))^2

@inline inv_ordinate(t::T, p::Int=1) where {T<:Real} = (asinh(log((one(T) + t) / (one(T) - t)) / T(π)))^(one(T) / T(p))

# Find the step size `h` at which a single step suffices to exhaust
# the numerical precision. Look at https://arxiv.org/pdf/2007.15057.pdf
find_tmax(T::Type, p::Int) = find_tmaxND(1, T, p)

# w == 0
# abs(x) == 1
# is that good?
function Level(h::T, n) where {T<:Real}
    points = Point{T}[]
    for i in 0:(2^(n - 1) - 1)
        t = h + i * (2h)
        x = ordinate(t)
        w = weight(t)
        # i > 0 && abs(x) == 1 && break
        # i > 0 && w == 0 && break
        push!(points, Point{T}(x, w))
    end
    reverse!(points)            # improve summation accuracy
    return Level{T}(points)
end

function TSQuadrature{T}(nlevels::Int=20, p::Int=1) where {T<:Real}
    @assert isodd(p)
    h = find_tmax(T, p)
    levels = Level{T}[]
    for level in 1:nlevels
        push!(levels, Level(h / 2^level, level))
    end
    return TSQuadrature{T}(h, levels)
end

# maps [a,b] to [-1,1]
transform(u::T, a::T, b::T) where {T} = (b + a) / T(2) + ((b - a) / T(2)) * u

function quadts(f, quad::TSQuadrature{T}, xmin::T, xmax::T; atol::T=zero(T),
                rtol::T=atol > 0 ? zero(T) : sqrt(eps(T))) where {T<:Real}
    if xmin == xmax
        return 0
    end
    Δx = (xmax - xmin) / 2
    h = quad.h

    x = (xmin + xmax) / 2
    w = weight(zero(T))
    s = h * w * f(transform(zero(T), xmin, xmax))
    levels = 0
    error = T(Inf) # norm(typeof(s)(Inf))

    for level in quad.levels
        h /= 2
        levels += 1
        sold = s
        s = zero(sold)
        for p in level.points
            xm = transform(-p.x, xmin, xmax)
            xp = transform(p.x, xmin, xmax)
            w = p.w
            # @show xm, xp, f(xm), f(xp), p.x
            s += w * (f(xm) + f(xp))
        end
        s = h * s + sold / 2

        tol = max(norm(s) * rtol, atol)
        error = norm(s - sold)
        levels ≥ 1 && error ≤ tol && break
    end

    return (result=Δx * s, error=error, levels=levels)
end

function quadts(f, quad::TSQuadrature{T}, xmin::Real, xmax::Real; atol::Real=zero(T),
                rtol::Real=atol > 0 ? zero(T) : sqrt(eps(T))) where {T<:Real}
    return quadts(f, quad, T(xmin), T(xmax); atol=T(atol), rtol=T(rtol))
end

################################################################################

export TSQuadratureND

# An integration point, defined by its coordinate `x` and weight `w`
struct PointND{D,T}
    x::SVector{D,T}
    w::T
end

# We integrate in levels. Each successivel level refines the previous
# one, and only contains the additional points and weights.
struct LevelND{D,T}
    points::Vector{PointND{D,T}}
end

# A complete quadrature scheme. `h` is the base step size for level 0.
struct TSQuadratureND{D,T}
    h::T
    levels::Vector{LevelND{D,T}}
end

# is w ≠ 0 ok?
function LevelND{D,T}(level::Level{T}) where {D,T<:Real}
    D::Int
    @assert D ≥ 0
    points = PointND{D,T}[]
    npoints = length(level.points)
    for i in CartesianIndex(ntuple(d -> 1, D)):CartesianIndex(ntuple(d -> npoints, D))
        x = SVector{D,T}(level.points[Tuple(i)[d]].x for d in 1:D)
        w = prod(level.points[Tuple(i)[d]].w for d in 1:D)
        # w ≠ 0 && push!(points, PointND{D,T}(x, w))
        push!(points, PointND{D,T}(x, w))
    end
    return LevelND{D,T}(points)
end

function find_tmaxND(D::Int, T::Type, p::Int)
    Fmin = eps(T)
    tmaxx = inv_ordinate(one(T) - Fmin, p)

    f(x, l) = weight(x, p)^D - Fmin
    u0 = (tmaxx)^(1 / D)
    probN = NonlinearProblem(f, u0)
    tmaxw = solve(probN, SimpleNewtonRaphson(); abstol=eps(T))[1]
    tmax = min(tmaxx, tmaxw)
    @show tmaxx, tmaxw
    return tmax
end

#if I add type in constructor, julia doesn't see the whole function, why?
function TSQuadratureND{D,T}(nlevels::Int, p::Int=1) where {D,T<:Real}
    D::Int #why ?
    @assert D ≥ 0

    h = find_tmaxND(D, T, p)
    levels = Level{T}[]
    for level in 1:nlevels
        push!(levels, Level(h / 2^level, level))
    end
    quad = TSQuadrature{T}(h, levels)
    levels = LevelND{D,T}[]
    for level in 1:nlevels
        push!(levels, LevelND{D,T}(quad.levels[level]))
        npoints = length(levels[end].points)
    end
    return TSQuadratureND{D,T}(h, levels)
end

# TSQuadratureND{D,T}(nlevels::Int) where {D,T<:Real} = TSQuadratureND{D,T}(TSQuadrature{T}(nlevels))
import Base.eltype
eltype(p::PointND{D,T}) where {D,T} = T
eltype(q::TSQuadratureND{D,T}) where {D,T} = T
eltype(l::LevelND{D,T}) where {D,T} = T

# maps [a,b] to [-1,1] in 3D
transform(u::AbstractVector{T}, a::AbstractVector{T}, b::AbstractVector{T}) where {T} = @. (b + a) / T(2) + ((b - a) / T(2)) * u

function quadts(f, quad::TSQuadratureND{D,T}, xmin::SVector{D,T}, xmax::SVector{D,T}; atol::T=zero(T),
                rtol::T=atol > 0 ? zero(T) : sqrt(eps(T))) where {D,T<:Real}
    D::Int
    @assert D ≥ 0

    Δx = (xmax - xmin) / 2
    h = quad.h # * Δx
    x = (xmin + xmax) / 2
    w = weight(zero(T))^D
    s = prod(h) * w * f(x...)

    levels = 0
    error = T(Inf) #norm(typeof(s)(Inf))
    for level in quad.levels
        h /= 2
        levels += 1
        sold = s

        s = zero(sold)
        for p in level.points
            t = zero(sold)
            for i in CartesianIndex(ntuple(d -> -1, D)):CartesianIndex(ntuple(d -> 2, D)):CartesianIndex(ntuple(d -> +1, D))
                xm = transform(-p.x, xmin, xmax)
                xp = transform(p.x, xmin, xmax)
                x = SVector{D,T}(Tuple(i)[d] < 0 ? xm[d] : xp[d] for d in 1:D)
                t += f(x...)
            end
            # p.w = w1*w2*w3...
            s += p.w * t
        end
        # shouldn't sold/2^D?
        s = h^D * s + sold / 2

        tol = max(norm(s) * rtol, atol)
        error = norm(s - sold)
        levels ≥ 1 && error ≤ tol && break
    end

    return (result=prod(Δx) * s, error=error, levels=levels)
end

export myquad2, myquad3
function myquad3(f, quad::TSQuadrature{T}, xmin::AbstractVector{T}, xmax::AbstractVector{T}) where {T<:Real}
    f1(x, y) = quadts(z -> f(x, y, z), quad, xmin[3], xmax[3])[1]
    f2(x) = quadts(y -> f1(x, y), quad, xmin[2], xmax[2])[1]
    f3() = quadts(x -> f2(x), quad, xmin[1], xmax[1])[1]
    return f3()
end

function myquad2(f, quad, xmin, xmax)
    f2(x) = quadts(y -> f(x, y), quad, xmin[2], xmax[2])[1]
    f3() = quadts(x -> f2(x), quad, xmin[1], xmax[1])[1]
    return f3()
end

function quadts(f, quad::TSQuadratureND{D,T}, xmin::AbstractVector{<:Real}, xmax::AbstractVector{<:Real};
                atol::Real=zero(T), rtol::Real=atol > 0 ? zero(T) : sqrt(eps(T))) where {D,T<:Real}
    return quadts(f, quad, SVector{D,T}(xmin), SVector{D,T}(xmax); atol=atol, rtol=rtol)
end

function quadts(f, quad::TSQuadratureND{D,T}, xmin::Union{NTuple{D},SVector{D}}, xmax::Union{NTuple{D},SVector{D}};
                atol::Real=zero(T), rtol::Real=atol > 0 ? zero(T) : sqrt(eps(T))) where {D,T<:Real}
    return quadts(f, quad, SVector{D,T}(xmin), SVector{D,T}(xmax); atol=atol, rtol=rtol)
end

end
