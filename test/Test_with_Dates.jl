using SchumakerSpline
using Dates
tol = 10*eps()

StartDate = Date(2018, 7, 21)
x = Array{Date}(1000)
for i in 1:1000
    x[i] = StartDate +Dates.Day(2* (i-1))
end

function f(x::Date)
    days_between = Dates.days(x - StartDate)
    return log(days_between+1) + sqrt(days_between)
end
y = f.(x)

spline = Schumaker(x,y)
for i in 1:length(x)
    abs(evaluate(spline, x[i]) - y[2]) < tol
end
# Evaluation with a float.
evaluate(spline, 11.5)

# Testing second derivatives
second_derivatives = evaluate.(spline, x,2)
maximum(second_derivatives) < tol

# Testing Integrals
function analytic_integral(lhs,rhs)
    lhs_in_days = Dates.days(lhs - StartDate)
    rhs_in_days = Dates.days(rhs - StartDate)
    return (rhs_in_days+1)*log(rhs_in_days+1)-rhs_in_days + (2/3)*rhs_in_days^(3/2) - ((lhs_in_days+1)*log(lhs_in_days+1) - lhs_in_days + (2/3)*lhs_in_days^(3/2))
end

lhs = StartDate
rhs = StartDate + Dates.Month(16)
numerical_integral = evaluate_integral(spline, lhs,rhs)
analytical = analytic_integral(lhs,rhs)
abs(  analytical - numerical_integral  ) < 1
