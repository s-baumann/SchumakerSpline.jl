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
########################
# Linear Extrapolation #
spline = Schumaker(x,y; extrapolation = (Linear, Linear))
# Now plotting the spline with derivatives
plot(spline, (-5.0, 10.0); derivs = true)
```

As a convenience the evaluate function can also be called with the shorthand:

```
vals = spline.(xrange)
derivative_values = spline.(xrange, 1)
second_derivative_values = spline.(xrange , 2)
```

We can now do the same with constant extrapolation.

```
##########################
# Constant Extrapolation #
spline = Schumaker(x,y; extrapolation = (Constant, Constant))
# Now plotting the spline with derivatives
plot(spline, (-5.0, 10.0); derivs = true)
```


If we did have gradient information we could get a better approximation by using it. In this case our gradients are:
```
analytical_first_derivative(e) = 1/e + 0.5 * e^(-0.5)
first_derivs = analytical_first_derivative.(x)
```
and we can generate a spline using these gradients with:
```
spline = Schumaker(x,y; gradients = first_derivs)
```
We could also have only specified the left or the right gradients using the left\_gradient and right\_gradient optional arguments.
