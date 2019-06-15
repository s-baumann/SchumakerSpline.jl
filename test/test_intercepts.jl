using SchumakerSpline
from = 0.5
to = 10
x1 = collect(range(from, stop=to, length=40))
y1 = (x1).^2
x2 = [0.5,0.75,0.8,0.93,0.9755,1.0,1.1,1.4,2.0]
y2 = sqrt.(x2)
s1 = Schumaker(x1,y1)
s2 = Schumaker(x2,y2)

crossover_point = get_intersection_points(s1,s2)
abs(evaluate(s1, crossover_point[1]) - evaluate(s2, crossover_point[1]))  < 100*eps()

y1 = (x1 .- 5).^2
x2 = [0.5,0.75,0.8,0.93,0.9755,1.0,1.1,1.4,2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
y2 = 1 .+ 0.5 .* x2
s1 = Schumaker(x1,y1)
s2 = Schumaker(x2,y2)

crossover_points = get_intersection_points(s1,s2)
abs(evaluate(s1, crossover_points[1]) - evaluate(s2, crossover_points[1]))  < 100*eps()
abs(evaluate(s1, crossover_points[2]) - evaluate(s2, crossover_points[2]))  < 100*eps()


y1 = (x1 .- 5).^2
x2 = [0.5,0.75,0.8,0.93,0.9755,1.0,1.1,1.4,2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
y2 = -0.5 .* x2
s1 = Schumaker(x1,y1)
s2 = Schumaker(x2,y2)

crossover_points = get_intersection_points(s1,s2)
length(crossover_points) == 0
