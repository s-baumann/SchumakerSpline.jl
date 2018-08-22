x = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]
y = log.(x) + sqrt.(x)
gradients = "Not-Supplied"

using SchumakerSpline: Schumaker
using Gadfly
using DataFrames
########################
# Linear Extrapolation #
extrapolation = "Linear"
spline = Schumaker(x,y; extrapolation = extrapolation)

df = DataFrame()
df[:x] =  range(-5, stop=10, length=100)
df[:values]  = evaluate.(spline, df[:x]  )
plt1 = Gadfly.plot(df, x=:x, y=:values, Geom.line)
df[:derivative_values]  = evaluate.(spline, df[:x] , 1 )
plt2 = Gadfly.plot(df, x=:x, y=:derivative_values, Geom.line)
df[:second_derivative_values]  = evaluate.(spline, df[:x] , 2 )
plt3 = Gadfly.plot(df, x=:x, y=:second_derivative_values, Geom.line)


##########################
# Constant Extrapolation #
extrapolation = "Constant"
spline = Schumaker(x,y; extrapolation = extrapolation)

df = DataFrame()
df[:x] =  range(-5, stop=10, length=100)
df[:values]  = evaluate.(spline, df[:x]  )
plt1 = Gadfly.plot(df, x=:x, y=:values, Geom.line)
df[:derivative_values]  = evaluate.(spline, df[:x] , 1 )
plt2 = Gadfly.plot(df, x=:x, y=:derivative_values, Geom.line)
df[:second_derivative_values]  = evaluate.(spline, df[:x] , 2 )
plt3 = Gadfly.plot(df, x=:x, y=:second_derivative_values, Geom.line)
##########################
# Curve    Extrapolation #
spline = Schumaker(x,y)

df = DataFrame()
df[:x] =  range(-5, stop=10, length=100)
df[:values]  = evaluate.(spline, df[:x]  )
plt1 = Gadfly.plot(df, x=:x, y=:values, Geom.line)
df[:derivative_values]  = evaluate.(spline, df[:x] , 1 )
plt2 = Gadfly.plot(df, x=:x, y=:derivative_values, Geom.line)
df[:second_derivative_values]  = evaluate.(spline, df[:x] , 2 )
plt3 = Gadfly.plot(df, x=:x, y=:second_derivative_values, Geom.line)
