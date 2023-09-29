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




## Related work

One-dimensional tanh-sinh quadrature:
[DoubleExponentialFormulas.jl](https://github.com/machakann/DoubleExponentialFormulas.jl)
