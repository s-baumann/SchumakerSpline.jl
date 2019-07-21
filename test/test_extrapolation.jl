using SchumakerSpline
from = 0.5
to = 10
x1 = collect(range(from, stop=to, length=40))
y1 = (x1).^2
x2 = [0.5,0.75,0.8,0.93,0.9755,1.0,1.1,1.4,2.0]
y2 = sqrt.(x2)
s1 = Schumaker(x1,y1; extrapolation = (Constant, Curve))
#plot(s1)
abs(s1(x1[1]) - s1(x1[1]-0.5)) < eps()
abs(s1(x1[40]) - s1(x1[40]+0.5)) > 0.01

s2 = Schumaker(x1,y1; extrapolation = (Constant, Constant))
abs(s2(x1[40]) - s2(x1[40]+0.5)) < eps()
plot(s2)

s3 = Schumaker(x1,y1; extrapolation = (Constant, Linear))
abs((s3(x1[40]+0.5) - s3(x1[40])) - (s3(x1[40]+1.5) - s3(x1[40]+1.0))) < 2e-14
