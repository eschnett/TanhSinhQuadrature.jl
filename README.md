# TanhSinhQuadrature.jl

Multi-dimensional [tanh-sinh
quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature) in
Julia

* [Documentation](https://eschnett.github.io/TanhSinhQuadrature.jl/dev/):
  Future documentation
* [GitHub](https://github.com/eschnett/TanhSinhQuadrature.jl): Source
  code repository
* [![GitHub
  CI](https://github.com/eschnett/TanhSinhQuadrature.jl/workflows/CI/badge.svg)](https://github.com/eschnett/TanhSinhQuadrature.jl/actions)
* [![codecov](https://codecov.io/gh/eschnett/TanhSinhQuadrature.jl/branch/main/graph/badge.svg?token=vHtLZhZpKG)](https://codecov.io/gh/eschnett/TanhSinhQuadrature.jl)

THIS PACKAGE IS INCOMPLETE. DO NOT USE.

I added conditions for the integration range $$[-t_n, t_n]$$ from this [paper](https://arxiv.org/pdf/2007.15057.pdf)
it gives improved error rates

``` julia
julia> f(x::T) where T = sqrt(T(2))/sqrt(one(T)+x)

julia> quadts(f, TSQuadrature{Float64}(), -1,1)
(result = 3.9999999347276813, error = 2.6538651720642292e-8, levels = 6)

julia> quadts(f, TSQuadrature{Float64}(), -1,1, atol=1e-16)
(result = 3.9999999782382543, error = 2.3940849303016876e-12, levels = 20)
```

Also the transformation now supports higher order of odd powers of the form 
$$\Psi(t) = \tanh(\frac{\pi}{2} \sinh(t^p)) $$
These are called TSHi where i is 1,3,5... 
for more details look into this [paper](https://ems.press/content/serial-article-files/2719)
``` julia
julia> f(x::T) where T = sqrt(T(2))/sqrt(one(T)+x)
# 1 is here corresponds to p=1 the usual TSH rule
julia> quadts(f, TSQuadrature{Float64}(1), -1,1)
(result = 3.9999999347276813, error = 2.6538651720642292e-8, levels = 6)

# 3 here is p=3 (TSH3 rule)
julia> quadts(f, TSQuadrature{Float64}(3), -1,1, atol=1e-16)
(result = 3.9999991492059976, error = 8.290316966252931e-7, levels = 20)
```


## issues

when the function has a singularity at one of the interval points,
and the interval is small, Underflow happens. As the levels of refinement
increase, the integration points coincide on t=+-1 and thus evaluate the
function at the singularity. 

This underflow instability stems from the variable transformation to map an
arbitrary domain [a,b] to [-1,1]. In the original conditions to avoid underflow 
instabilities from this [paper](https://arxiv.org/pdf/2007.15057.pdf), the authors
don't take into account this variable trasformation. We provide a function 
to find $t_max$ taking into account the variable transformation, but the need to specify
the domain of integration defeats the design of a reusable {TSQuadrature} object.

A workaround is to use higher order TSHi methods. For example TSH3 works and avoids
the underflow instability due to the more rapid fall off of the transformed function.
``` julia
julia> f(x::T) where T = one(T)/(x-one(T))

# 1 here corresponds to p=1, the original TSH rule
julia> quad = TSQuadrature{T}(10, 1);

julia> quadts(h, quad, 0,1)
(-36.64917027223108, 0.09679325694600749, 10)

julia> quadts(h, quad, 0.92,1)
(-34.66297119486927, 0.01782170555074116, 10)

julia> quadts(h, quad, 0.93,1)
(Inf, Inf, 6)

# 3 here corresponds to p=3, the TSH3 rule
julia> quad = TSQuadrature{Float64}(10, 3);

julia> quadts(h, quad, 0,1)
(-6.440474048946902, 0.005131276307238863, 10)

julia> quadts(h, quad, 0.92,1)
(-6.44047404894647, 0.005131276307250798, 10)

julia> quadts(h, quad, 0.93,1)
(-6.440474048947368, 0.0051312763072301915, 10)

julia> quadts(h, quad, 0.999999999,1)
(-6.440440508722164, 0.0051306262629217175, 10)
```

TL;DR : if you get Inf use p>=3 (odd)

## Related work

One-dimensional tanh-sinh quadrature:
[DoubleExponentialFormulas.jl](https://github.com/machakann/DoubleExponentialFormulas.jl)
