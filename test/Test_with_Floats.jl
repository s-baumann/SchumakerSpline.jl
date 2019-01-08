using SchumakerSpline
tol = 10*eps()

x = collect(range(1, stop=6, length=1000))
y = log.(x) + sqrt.(x)

spline = Schumaker(x,y)
for i in 1:length(x)
    abs(evaluate(spline, x[i]) - y[2]) < tol
end

# Testing First derivative.
first_derivatives = evaluate.(spline, x, 1)
analytical_first_derivative(e) = 1/e + 0.5 * e^(-0.5)
maximum(abs.(first_derivatives .- analytical_first_derivative.(x))) < 0.002

# Testing second derivatives
second_derivatives = evaluate.(spline, x, 2)
maximum(second_derivatives) < tol

# Test higher derivative
higher_derivatives = evaluate.(spline, x, 3)
maximum(second_derivatives) < tol

# Testing Integrals
analytic_integral(lhs,rhs) = rhs*log(rhs) - rhs + (2/3) * rhs^(3/2) - ( lhs*log(lhs) - lhs + (2/3) * lhs^(3/2) )
lhs = 2.0
rhs = 2.5
numerical_integral = evaluate_integral(spline, lhs,rhs)
abs(analytic_integral(lhs,rhs) - numerical_integral) < 0.01

lhs = 1.2
rhs = 4.3
numerical_integral = evaluate_integral(spline, lhs,rhs)
abs(analytic_integral(lhs,rhs) - numerical_integral) < 0.01

lhs = 0.8
rhs = 4.0
numerical_integral = evaluate_integral(spline, lhs,rhs)
abs(analytic_integral(lhs,rhs) - numerical_integral) < 0.03
