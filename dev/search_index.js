var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "FixedPointAcceleration.jl",
    "title": "FixedPointAcceleration.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#FixedPointAcceleration.jl-1",
    "page": "FixedPointAcceleration.jl",
    "title": "FixedPointAcceleration.jl",
    "category": "section",
    "text": "Fixed point finders are conceptually similar to both optimisation and root finding algorithms but thus far implementations of fixed point finders have been rare in Julia. In some part this is likely because there is often an obvious method to find a fixed point by merely feeding a guessed fixed point into a function, taking the result and feeding it back into the function. By doing this repeatedly a fixed point is often found. This method (that we will call the \"Simple\" method) is often convergent but it is also often slow which can be prohibitive when the function itself is expensive.FixedPointAcceleration.jl aims to provide fixed point acceleration algorithms that can be much faster than the simple method. In total eight algorithms are implemented. The first is the simple method as described earlier. There are also the Newton, Aitken and Scalar Epsilon Algorithm (SEA) methods that are designed for accelerating the convergence of scalar sequences. Four algorithms for accelerating vector sequences are also implemented including the Vector Epsilon Algorithm (VEA), two minimal polynomial algorithms (MPE and RRE)  and Anderson acceleration.In this paper section 1 starts by with a brief explanation of fixed points before section 2 describes the acceleration algorithms provided by FixedPointAcceleration.jl. Here the goal is  to illustrate the logic underling each algorithm so users can better choose which suits their problem. Readers interested in the underlying mathematics are referred to the original papers. Section 3 then illustrates a few features of the package that enable users to better track the progress of an algorithm while it is running and switch algorithms if desired before a fixed point is found.Section 4 then presents several applications of these fixed point algorithms in economics, asset pricing and machine learning. Finally section 5 presents a convergence comparison showing how many iterations each algorithm takes in bringing each problem to its fixed point for each of the applications presented in section 4.pages = [\"index.md\",\n         \"examples.md\"]\nDepth = 2"
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
