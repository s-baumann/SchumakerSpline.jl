using Documenter, SchumakerSpline

makedocs(modules=[SchumakerSpline],
)

deploydocs(deps = Deps.pip("pygments", "mkdocs", "mkdocs-material", "python-markdown-math"),
    repo = "github.com/s-baumann/SchumakerSpline.jl.git",
    julia = "0.4",
    osname = "windows",
)
