using Documenter, SchumakerSpline

makedocs(
    format = Documenter.HTML(),
    sitename = "SchumakerSpline",
    modules = [SchumakerSpline],
    pages = Any[
        "Introduction" => "index.md",
        "Examples" => "examples.md",
        "API" => "api.md"]
    ]
)

deploydocs(
    repo   = "github.com/s-baumann/SchumakerSpline.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
