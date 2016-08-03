using SchumakerSpline

x = [1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]
y = log(x) + sqrt(x)

Spline, DerivSpline, SecondDerivSpline = schumaker(x,y)

xArray = linspace(1, 6, 100)
Spp = SecondDerivSpline(xArray)

maximum(Spp) < 2*eps()
