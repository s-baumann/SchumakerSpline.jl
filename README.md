# SchumakerSpline

[![Build Status](https://travis-ci.org/s-baumann/SchumakerSpline.jl.svg?branch=master)](https://travis-ci.org/s-baumann/SchumakerSpline.jl)

A Julia package to create a shape preserving spline. This is guaranteed to be monotonic and concave or convex if the data is monotonic and concave or convex. It does not use any optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. This package has the same functionality as the R package called [schumaker](https://cran.r-project.org/web/packages/schumaker/index.html).

See the docs [here](https://s-baumann.github.io/SchumakerSpline.jl/)

Note that this package currently supports Julia 1.0 and 0.7 but potentially not earlier versions of Julia. Authored by Stuart Baumann and Margaryta Klymak this software is available to all under the standard MIT license. See LICENSE.md
