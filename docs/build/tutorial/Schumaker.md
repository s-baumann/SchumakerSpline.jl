
<a id='schumaker-Function-1'></a>

# schumaker Function


The only exported function in this package (and the only one a basic user will ever need) is the schumaker function which creates a spline along with functions giving its derivative and second derivative. This is described as below


<a id='Function-documentation-1'></a>

## Function documentation

<a id='SchumakerSpline.schumaker' href='#SchumakerSpline.schumaker'>#</a>
**`SchumakerSpline.schumaker`** &mdash; *Function*.



Creates splines for a given set of x and y values (and optionally gradients) and the first and second derivatives of this spline.

**Takes**

  * x - A Float64 vector of x coordinates.
  * y - A Float64 vector of y coordinates.
  * gradients (optional)- A Float64 vector of gradients at each point. If not supplied these are imputed from x and y.
  * extrapolation (optional) - This should be a string in ("Curve", "Linear", "Constant") specifying how to interpolate outside of the sample domain. By default it is "curve" which extends out the first and last quadratic curves. The other options are "Linear" which extends the line (from first and last curve) out from the first and last point and "Constant" which extends out the y value at the first and last point. ### Returns
  * A spline which takes and input value and returns the spline y value.
  * The derivative of this spline.
  * The second derivative of this spline.  

