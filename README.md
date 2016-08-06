# SchumakerSpline

Linux, OSX: [![Build Status](https://travis-ci.org/s-baumann/SchumakerSpline.jl.svg?branch=master)](https://travis-ci.org/s-baumann/SchumakerSpline.jl)
Windows: [![Build status](https://ci.appveyor.com/api/projects/status/o7g0etv2i5udghj8?svg=true)](https://ci.appveyor.com/project/s-baumann/schumakerspline-jl)

A Julia package to create a shape preserving spline. This is guaranteed to be monotonic and concave or convex if the data is monotonic and concave or convex. It does not use any optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. It also automatically gives the first two derivatives
of the spline and options for determining behaviour when evaluated outside the interpolation domain.

This package has the same functionality as the R package called [schumaker](https://cran.r-project.org/web/packages/schumaker/index.html).

Note that this package does not currently support Julia 0.5 but will be amended to support it after it becomes the stable version. Authored by Stuart Baumann and Margaryta Klymak this software is available to all under the standard MIT license. See LICENSE.md
