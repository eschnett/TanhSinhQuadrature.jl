using DoubleFloats
using Random
using StaticArrays
using TanhSinhQuadrature
using Test

################################################################################

const Types = [Float32, Float64, Double64, BigFloat]

################################################################################

const quads = Dict()

@testset "Create quadrature rules T=$T" for T in Types
    quads[T] = TSQuadrature{T}()
end

@testset "Basic integration T=$T" for T in Types
    quad = quads[T]::TSQuadrature{T}

    f0(x) = 1
    f1(x) = x + 1
    f2(x) = 3 * x^2

    @test quadts(f0, quad, -1, +1).result ≈ 2
    @test quadts(f1, quad, -1, +1).result ≈ 2
    @test quadts(f2, quad, -1, +1).result ≈ 2
end

Random.seed!(0)
@testset "Integral bounds T=$T" for T in Types
    quad = quads[T]::TSQuadrature{T}

    f0(x) = 1
    f1(x) = x + 1
    f2(x) = 3 * x^2

    xmin = -1 + T(rand(-5:5)) / 10
    xmax = +1 + T(rand(-5:5)) / 10

    @test quadts(f0, quad, xmin, xmax).result ≈ xmax - xmin
    @test quadts(f1, quad, xmin, xmax).result ≈ (xmax^2 - xmin^2) / 2 + xmax - xmin
    @test quadts(f2, quad, xmin, xmax).result ≈ (xmax^3 - xmin^3)
end

Random.seed!(0)
@testset "Linearity T=$T" for T in Types
    quad = quads[T]::TSQuadrature{T}

    a = 1 + T(rand(-5:5)) / 10
    b0 = 1 + T(rand(-5:5)) / 10
    b1 = 1 + T(rand(-5:5)) / 10
    b2 = 1 + T(rand(-5:5)) / 10
    c0 = 1 + T(rand(-5:5)) / 10
    c1 = 1 + T(rand(-5:5)) / 10
    c2 = 1 + T(rand(-5:5)) / 10

    f(x) = b0 + b1 * x + b2 * x^2
    g(x) = c0 + c1 * x + c2 * x^2

    F = quadts(f, quad, -1, +1).result
    G = quadts(g, quad, -1, +1).result

    afg(x) = a * f(x) + g(x)
    @test quadts(afg, quad, -1, +1).result ≈ a * F + G

    d = T(rand(-9:9)) / 10

    @test quadts(f, quad, +1, -1).result ≈ -F

    F0 = quadts(f, quad, -1, a).result
    F1 = quadts(f, quad, a, +1).result
    @test F0 + F1 ≈ F
end

Random.seed!(0)
@testset "Integrals of singular functions T=$T" for T in Types
    quad = quads[T]::TSQuadrature{T}

    rtol = Dict(Float32 => 10 * sqrt(eps(T)),
                Float64 => 10 * sqrt(eps(T)),
                Double64 => 10^5 * sqrt(eps(T)),
                BigFloat => 1.0e-21)[T]

    f(x) = 1 / sqrt(1 - x^2)
    F = quadts(f, quad, -1, +1; rtol=eps(T)^(T(3) / 4)).result
    @test isapprox(F, T(π); rtol=rtol)

    rtol = Dict(Float32 => 10 * sqrt(eps(T)),
                Float64 => 100 * sqrt(eps(T)),
                Double64 => 10^7 * sqrt(eps(T)),
                BigFloat => 1.0e-18)[T]

    a = T(rand(1:9)) / 10
    g(x) = 1 / x^(1 - a)
    G = quadts(g, quad, 0, 1; rtol=eps(T)^(T(3) / 4)).result
    @test isapprox(G, 1 / a; rtol=rtol)
end
