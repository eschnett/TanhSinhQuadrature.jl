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
function find_tmax(T::Type, p::Int)
    Fmin = eps(T)
    tmaxx = inv_ordinate(one(T) - Fmin, p)

    # y(x, l) = ordinate(x, p) - 1 + Fmin
    # u0 = tmaxx
    # probN = NonlinearProblem(y, u0)
    # tmaxx = solve(probN, SimpleNewtonRaphson(); abstol=1e-20)[1]

    f(x, l) = weight(x, p) - Fmin
    u0 = tmaxx
    probN = NonlinearProblem(f, u0)
    tmaxw = solve(probN, SimpleNewtonRaphson(); abstol=eps(T))[1]
    @show tmaxw
    @show tmaxx
    return min(tmaxx, tmaxw)
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

function TSQuadrature{T}(p::Int=1, nlevels::Int=20) where {T<:Real}
    @assert isodd(p)
    h = find_tmax(T, p)
    @show h
    levels = Level{T}[]
    # println("TSQuadrature{$T}:")
    for level in 1:nlevels
        push!(levels, Level(h / 2^level))
        npoints = length(levels[end].points)
        # println("    level: $level    npoints: $npoints")
    end
    return TSQuadrature{T}(h, levels)
end

function quadts(f, quad::TSQuadrature{T}, xmin::T, xmax::T; atol::T=zero(T),
                rtol::T=atol > 0 ? zero(T) : sqrt(eps(T))) where {T<:Real}
    Δx = (xmax - xmin) / 2
    h = quad.h * Δx

    x = (xmin + xmax) / 2
    w = weight(zero(T))
    s = h * w * f(x)
    levels = 0
    error = norm(typeof(s)(Inf))

    for level in quad.levels
        h /= 2
        levels += 1
        sold = s

        s = zero(sold)
        for p in level.points
            xm = xmin + Δx * (1 - p.x)
            xp = xmax - Δx * (1 - p.x)
            w = p.w
            s += w * (f(xm) + f(xp))
        end
        s = h * s + sold / 2

        tol = max(norm(s) * rtol, atol)
        error = norm(s - sold)
        levels ≥ 1 && error ≤ tol && break
    end

    return (result=s, error=error, levels=levels)
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

function LevelND{D,T}(level::Level{T}) where {D,T<:Real}
    D::Int
    @assert D ≥ 0
    points = PointND{D,T}[]
    npoints = length(level.points)
    for i in CartesianIndex(ntuple(d -> 1, D)):CartesianIndex(ntuple(d -> npoints, D))
        x = SVector{D,T}(level.points[Tuple(i)[d]].x for d in 1:D)
        w = prod(level.points[Tuple(i)[d]].w for d in 1:D)
        w ≠ 0 && push!(points, PointND{D,T}(x, w))
    end
    return LevelND{D,T}(points)
end

function TSQuadratureND{D,T}(quad::TSQuadrature{T}) where {D,T<:Real}
    D::Int
    @assert D ≥ 0
    h = quad.h
    nlevels = length(quad.levels)
    levels = LevelND{D,T}[]
    # println("TSQuadratureND{$D,$T}:")
    for level in 1:nlevels
        push!(levels, LevelND{D,T}(quad.levels[level]))
        npoints = length(levels[end].points)
        # println("    level: $level    npoints: $npoints")
    end
    return TSQuadratureND{D,T}(h, levels)
end

TSQuadratureND{D,T}(nlevels::Int) where {D,T<:Real} = TSQuadratureND{D,T}(TSQuadrature{T}(nlevels))

function quadts(f, quad::TSQuadratureND{D,T}, xmin::SVector{D,T}, xmax::SVector{D,T}; atol::T=zero(T),
                rtol::T=atol > 0 ? zero(T) : sqrt(eps(T))) where {D,T<:Real}
    D::Int
    @assert D ≥ 0

    Δx = (xmax - xmin) / 2
    h = quad.h^D * prod(Δx)

    x = (xmin + xmax) / 2
    w = weight(zero(T))^D
    s = h * w * f(x...)
    levels = 0
    error = norm(typeof(s)(Inf))

    for level in quad.levels
        h /= 2^D
        levels += 1
        sold = s

        s = zero(sold)
        for p in level.points
            t = zero(s)
            for i in CartesianIndex(ntuple(d -> -1, D)):CartesianIndex(ntuple(d -> 2, D)):CartesianIndex(ntuple(d -> +1, D))
                xm = xmin + Δx .* (1 .- p.x)
                xp = xmax - Δx .* (1 .- p.x)
                x = SVector{D,T}(Tuple(i)[d] < 0 ? xm[d] : xp[d] for d in 1:D)
                t += f(x...)
            end
            w = p.w
            s += w * t
        end
        s = h * s + sold / 2

        tol = max(norm(s) * rtol, atol)
        error = norm(s - sold)
        levels ≥ 4 && error ≤ tol && break
    end

    return (result=s, error=error, levels=levels)
end

function quadts(f, quad::TSQuadratureND{D,T}, xmin::Union{NTuple{D},SVector{D}}, xmax::Union{NTuple{D},SVector{D}};
                atol::Real=zero(T), rtol::Real=atol > 0 ? zero(T) : sqrt(eps(T))) where {D,T<:Real}
    return quadts(f, quad, SVector{D,T}(xmin), SVector{D,T}(xmax); atol=atol, rtol=rtol)
end

end
