var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "SchumakerSpline.jl",
    "title": "SchumakerSpline.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#SchumakerSpline.jl-1",
    "page": "SchumakerSpline.jl",
    "title": "SchumakerSpline.jl",
    "category": "section",
    "text": "A Julia package to create a shape preserving spline. This is guaranteed to be monotonic and concave or convex if the data is monotonic and concave or convex. It does not use any optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. This package has the same functionality as the R package called schumaker.pages = [\"index.md\",\n         \"examples.md\"]\nDepth = 2"
},

{
    "location": "examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "Generating some example datax = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]\ny = log.(x) + sqrt.(x)\ngradients = missingIn this case we do not have gradients information and so gradients will be imputed from the x and y data.We can create a spline and plot it with linear extrapolation.using SchumakerSpline\nusing Plots\n########################\n# Linear Extrapolation #\nspline = Schumaker(x,y; extrapolation = Linear)\n# Now plotting the spline\nxrange =  collect(range(-5, stop=10, length=100))\nvalues  = evaluate.(spline, xrange)\nderivative_values  = evaluate.(spline, xrange, 1 )\nsecond_derivative_values  = evaluate.(spline, xrange , 2 )\nplot(xrange , values; label = \"Spline\")\nplot!(xrange, derivative_values; label = \"First Derivative\")\nplot!(xrange, second_derivative_values; label = \"Second Derivative\")We can now do the same with constant extrapolation.##########################\n# Constant Extrapolation #\nextrapolation = Constant\nspline = Schumaker(x,y; extrapolation = Constant)\n# Now plotting the spline\nxrange =  collect(range(-5, stop=10, length=100))\nvalues  = evaluate.(spline, xrange)\nderivative_values  = evaluate.(spline, xrange, 1 )\nsecond_derivative_values  = evaluate.(spline, xrange , 2 )\nplot(xrange , values; label = \"Spline\")\nplot!(xrange, derivative_values; label = \"First Derivative\")\nplot!(xrange, second_derivative_values; label = \"Second Derivative\")"
},

]}
