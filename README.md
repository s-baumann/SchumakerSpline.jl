# SchumakerSpline

| Build | Coverage | Documentation |
|-------|----------|---------------|
| [![Build status](https://github.com/s-baumann/SchumakerSpline.jl/workflows/CI/badge.svg)](https://github.com/s-baumann/SchumakerSpline.jl/actions) | [![codecov](https://codecov.io/gh/s-baumann/SchumakerSpline.jl/branch/main/graph/badge.svg?token=sElLVJgRel)](https://codecov.io/gh/s-baumann/SchumakerSpline.jl) | [![docs-latest-img](https://img.shields.io/badge/docs-latest-blue.svg)](https://s-baumann.github.io/SchumakerSpline.jl/dev/index.html) |

A Julia package to create a shape preserving spline. This is guaranteed to be monotonic and concave or convex if the data is monotonic and concave or convex. It does not use any optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. This package has the same functionality as the R package called [schumaker](https://cran.r-project.org/web/packages/schumaker/index.html).
