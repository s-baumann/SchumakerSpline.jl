using Documenter
using SchumakerSpline

makedocs(
 # options
 modules=[SchumakerSpline],
 # html options
 format = :html,
 sitename = "SchumakerSpline.jl",
 pages = [
     "Introduction" => ["index.md"],
 ],
html_canonical = "https://s-baumann.github.io/SchumakerSpline.jl/latest/",
)
deploydocs(
 repo = "github.com/s-baumann/SchumakerSpline.jl.git",
 target = "build",
 julia = "1.0",
 deps = nothing,
 make = nothing
)
