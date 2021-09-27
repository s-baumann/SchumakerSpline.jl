using SchumakerSpline
using Test
tol = 100 * eps()
close(a,b) = abs(a-b) < tol

from = 0.5
to = 10
x1 = collect(range(from, stop=to, length=40))
y1 = (x1).^2
x2 = [0.5,0.75,0.8,0.93,0.9755,1.0,1.1,1.4,2.0]
y2 = sqrt.(x2)
s1 = Schumaker(x1,y1)
s2 = Schumaker(x2,y2)
con = 7.6
x = 0.9

close( (s1 +con)(x) , s1(x) + con)
close( (con + s1)(x) , s1(x) + con)

close( (s1 - con)(x) , s1(x) - con)
close( (con - s1)(x) , con - s1(x) )

close( (s1 *con)(x) , s1(x) * con)
close( (con * s1)(x) , s1(x) * con)

close( (s1 /con)(x) , s1(x) / con)
