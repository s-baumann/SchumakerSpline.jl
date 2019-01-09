# Examples

Generating some example data

```
x = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]
y = log.(x) + sqrt.(x)
gradients = missing
```
In this case we do not have gradients information and so gradients will be imputed from the x and y data.

We can create a spline and plot it with linear extrapolation.

```
using SchumakerSpline
using Plots
########################
# Linear Extrapolation #
spline = Schumaker(x,y; extrapolation = Linear)
# Now plotting the spline
xrange =  collect(range(-5, stop=10, length=100))
values  = evaluate.(spline, xrange)
derivative_values  = evaluate.(spline, xrange, 1 )
second_derivative_values  = evaluate.(spline, xrange , 2 )
plot(xrange , values; label = "Spline")
plot!(xrange, derivative_values; label = "First Derivative")
plot!(xrange, second_derivative_values; label = "Second Derivative")
```

We can now do the same with constant extrapolation.

```
##########################
# Constant Extrapolation #
extrapolation = Constant
spline = Schumaker(x,y; extrapolation = Constant)
# Now plotting the spline
xrange =  collect(range(-5, stop=10, length=100))
values  = evaluate.(spline, xrange)
derivative_values  = evaluate.(spline, xrange, 1 )
second_derivative_values  = evaluate.(spline, xrange , 2 )
plot(xrange , values; label = "Spline")
plot!(xrange, derivative_values; label = "First Derivative")
plot!(xrange, second_derivative_values; label = "Second Derivative")
```
