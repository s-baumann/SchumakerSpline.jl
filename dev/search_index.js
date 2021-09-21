var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = SchumakerSpline","category":"page"},{"location":"api/#Internal-Functions","page":"API","title":"Internal Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Pages = [\"api.md\"]","category":"page"},{"location":"api/#Main-Struct","page":"API","title":"Main Struct","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"    Schumaker\n    evaluate\n    evaluate_integral","category":"page"},{"location":"api/#SchumakerSpline.Schumaker","page":"API","title":"SchumakerSpline.Schumaker","text":"Schumaker(x::Array{T,1},y::Array{T,1} ; gradients::Union{Missing,Array{T,1}} = missing, extrapolation::Schumaker_ExtrapolationSchemes = Curve,\n              left_gradient::Union{Missing,T} = missing, right_gradient::Union{Missing,T} = missing)\nSchumaker(x::Array{Int,1},y::Array{T,1} ; gradients::Union{Missing,Array{T,1}} = missing, extrapolation::Schumaker_ExtrapolationSchemes = Curve,\n              left_gradient::Union{Missing,T} = missing, right_gradient::Union{Missing,T} = missing)\nSchumaker(x::Array{Date,1},y::Array{T,1} ; gradients::Union{Missing,Array{T,1}} = missing, extrapolation::Schumaker_ExtrapolationSchemes = Curve,\n              left_gradient::Union{Missing,T} = missing, right_gradient::Union{Missing,T} = missing)\n\nCreates a Schumaker spline.\n\nTakes\n\nx - A vector of x coordinates.\ny - A vector of y coordinates.\nextrapolation (optional) - This should be Curve, Linear or Constant specifying how to interpolate outside of the sample domain. By default it is Curve which extends out the first and last quadratic curves. The other options are Linear which extends the line (from first and last curve) out from the first and last point and Constant which extends out the y value at the first and last point.\ngradients (optional)- A vector of gradients at each point. If not supplied these are imputed from x and y.\nleft_gradient - The gradient at the lowest value of x in the domain. This will override the gradient imputed or submitted in the gradients optional argument (if it is submitted there)\nright_gradient - The gradient at the highest value of x in the domain. This will override the gradient imputed or submitted in the gradients optional argument (if it is submitted there)\n\nReturns\n\nA Schumaker object which details the spline. This object can then be evaluated with evaluate or evaluate_integral.\n\n\n\n\n\n\n\n","category":"type"},{"location":"api/#SchumakerSpline.evaluate","page":"API","title":"SchumakerSpline.evaluate","text":"Evaluates the spline at a point. The point can be specified as a Real number (Int, Float, etc) or a Date. Derivatives can also be taken.\n\nTakes\n\nspline - A Schumaker type spline\nPointToExamine - The point at which to evaluate the integral\nderivative - The derivative being sought. This should be 0 to just evaluate the spline, 1 for the first derivative or 2 for a second derivative.\n\nHigher derivatives are all zero (because it is a quadratic spline). Negative values do not give integrals. Use evaluate_integral instead.\n\nReturns\n\nA value of the spline or appropriate derivative in the same format as specified in the spline.\n\n\n\n\n\n","category":"function"},{"location":"api/#SchumakerSpline.evaluate_integral","page":"API","title":"SchumakerSpline.evaluate_integral","text":"Estimates the integral of the spline between lhs and rhs. These end points can be input as Reals or Dates.\n\nTakes\n\nspline - A Schumaker type spline\nlhs - The left hand limit of the integral\nrhs - The right hand limit of the integral\n\nReturns\n\nA Float64 value of the integral.\n\n\n\n\n\n","category":"function"},{"location":"api/#Extrapolation-Schemes","page":"API","title":"Extrapolation Schemes","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"    Schumaker_ExtrapolationSchemes\n    Curve\n    Linear\n    Constant","category":"page"},{"location":"api/#SchumakerSpline.Schumaker_ExtrapolationSchemes","page":"API","title":"SchumakerSpline.Schumaker_ExtrapolationSchemes","text":"This creates an enum which details how extrapolation from the interpolation domain should be done.\n\n\n\n\n\n","category":"type"},{"location":"api/#Working-with-Splines","page":"API","title":"Working with Splines","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"find_derivative_spline\nfind_roots\nfind_optima\nget_crossover_in_interval\nget_intersection_points\nreshape_values\nsplice_splines","category":"page"},{"location":"api/#SchumakerSpline.find_derivative_spline","page":"API","title":"SchumakerSpline.find_derivative_spline","text":"findderivativespline(spline::Schumaker)     Returns a SchumakerSpline that is the derivative of the input spline\n\n\n\n\n\n","category":"function"},{"location":"api/#SchumakerSpline.find_roots","page":"API","title":"SchumakerSpline.find_roots","text":"findroot(spline::Schumaker; rootvalue::T = 0.0)     Finds roots - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find roots.     Here root_value can be set to get all points at which the function is equal to the root value. For instance if you want to find all points at which     the spline has a value of 1.0.\n\n\n\n\n\n","category":"function"},{"location":"api/#SchumakerSpline.find_optima","page":"API","title":"SchumakerSpline.find_optima","text":"find_optima(spline::Schumaker) Finds optima - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find optima.\n\n\n\n\n\n","category":"function"},{"location":"api/#SchumakerSpline.get_crossover_in_interval","page":"API","title":"SchumakerSpline.get_crossover_in_interval","text":"get_crossover_in_interval(s1::Schumaker{T}, s2::Schumaker{R}, interval::Tuple{U,U}) where T<:Real where R<:Real where U<:Real\n\nFinds the point at which two schumaker splines cross over each other within a single interval. This is not exported.\n\n\n\n\n\n","category":"function"},{"location":"api/#Two-dimensional-splines","page":"API","title":"Two dimensional splines","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Schumaker2d","category":"page"},{"location":"api/#Plotting-of-splines","page":"API","title":"Plotting of splines","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"plot","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Generating some example data","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"x = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]\ny = log.(x) + sqrt.(x)\ngradients = missing","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"In this case we do not have gradients information and so gradients will be imputed from the x and y data.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can create a spline and plot it with linear extrapolation.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using SchumakerSpline\nusing Plots\n########################\n# Linear Extrapolation #\nspline = Schumaker(x,y; extrapolation = (Linear, Linear))\n# Now plotting the spline\nxrange = collect(range(-5, stop=10, length=100))\nvals = evaluate.(spline, xrange)\nderivative_values = evaluate.(spline, xrange, 1 )\nsecond_derivative_values = evaluate.(spline, xrange , 2 )\nplot(xrange , vals; label = \"Spline\")\nplot!(xrange, derivative_values; label = \"First Derivative\")\nplot!(xrange, second_derivative_values; label = \"Second Derivative\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"As a convenience the evaluate function can also be called with the shorthand:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"vals = spline.(xrange)\nderivative_values = spline.(xrange, 1)\nsecond_derivative_values = spline.(xrange , 2)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We can now do the same with constant extrapolation.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"##########################\n# Constant Extrapolation #\nextrapolation = (Constant, Constant)\nspline = Schumaker(x,y; extrapolation = extrapolation)\n# Now plotting the spline\nxrange =  collect(range(-5, stop=10, length=100))\nvals  = evaluate.(spline, xrange)\nderivative_values  = evaluate.(spline, xrange, 1 )\nsecond_derivative_values  = evaluate.(spline, xrange , 2 )\nplot(xrange , vals; label = \"Spline\")\nplot!(xrange, derivative_values; label = \"First Derivative\")\nplot!(xrange, second_derivative_values; label = \"Second Derivative\")","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"If we did have gradient information we could get a better approximation by using it. In this case our gradients are:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"analytical_first_derivative(e) = 1/e + 0.5 * e^(-0.5)\nfirst_derivs = analytical_first_derivative.(x)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"and we can generate a spline using these gradients with:","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"spline = Schumaker(x,y; gradients = first_derivs)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"We could also have only specified the left or the right gradients using the left_gradient and right_gradient optional arguments.","category":"page"},{"location":"#SchumakerSpline.jl","page":"Introduction","title":"SchumakerSpline.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"A simple shape preserving spline implementation in Julia.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"A Julia package to create a shape preserving spline. This is a shape preserving spline which is guaranteed to be monotonic and concave/convex if the data is monotonic and concave/convex. It does not use any numerical optimisation and is therefore quick and smoothly converges to a fixed point in economic dynamics problems including value function iteration. Analytical derivatives and integrals of the spline can easily be taken through the evaluate and evaluate_integral functions.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This package has the same basic functionality as the R package called schumaker.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"If you want to do algebraic operations on splines you can also use a schumaker spline through the UnivariateFunctions package.","category":"page"},{"location":"#Optional-parameters","page":"Introduction","title":"Optional parameters","text":"","category":"section"},{"location":"#Gradients.","page":"Introduction","title":"Gradients.","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"The gradients at each of the (x,y) points can be input to give more accuracy. If not supplied these are estimated from the points provided. It is also possible to input on the gradients on the edges of the x domain and have all of the intermediate gradients imputed.","category":"page"},{"location":"#Out-of-sample-prediction.","page":"Introduction","title":"Out of sample prediction.","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"There are three options for out of sample prediction.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Curve - This is where the quadratic curve that is present in the first and last interval are used to predict points before the first interval and after the last interval respectively.\nLinear - This is where a line is extended out before the first interval and after the last interval. The slope of the line is given by the derivative at the start of the first interval and end of the last interval.\nConstant - This is where the first and last y values are used for prediction before the first point of the interval and after the last part of the interval respectively.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"pages = [\"index.md\",\n         \"examples.md\"]\nDepth = 2","category":"page"}]
}
