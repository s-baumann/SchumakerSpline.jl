
<a id='schumaker-Function-1'></a>

# Schumaker Object

To use this package you first create a Schumaker object. The constructor for this takes a vector of x values, a vector of y values and optionally a vector of gradients and a specification of the type of extrapolation that is desired.
Then you can use the evaluate function or the evaluate_integral functions on this object. The evaluate function will evaluate the spline or one of the splines derivates at a point. The evaluate_integral function will give the integral of the spline between two points.

On the x domain, the constructor of the Schumaker object can take in either Float64s, Ints or Dates (from the Dates package). A FLoat64, Int or Date can similarly be input to the evaluate or evaluate_integral functions.

# Backwards compatibility

There is no backwards compatibility with the previous version of SchumakerSpline.jl. This package has been completely changed. It is now based around a type for which functions and types can be written (so it is now feasible to have the spline as a member object of a Forward_FX_Curve type). It now has support for dates as an x input. In addition it is now possible to integrate the spline analytically.
