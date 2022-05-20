# Run `julia --project=@. --color=yes make.jl` in this directory to
# build the documentation. The `Documenter` must be available (i.e.
# installed globally) for this to work.

using Documenter

using TanhSinhQuadrature

makedocs(; sitename="TanhSinhQuadrature", format=Documenter.HTML(; prettyurls=false))
deploydocs(; repo="github.com/eschnett/TanhSinhQuadrature.jl")
