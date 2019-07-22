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

# Testing Rootfinder and OptimaFinder
from = 0.0
to = 10.0
x = collect(range(from, stop=to, length=400))
# This should have no roots or optima.

y = log.(x) + sqrt.(x)
spline = Schumaker(x,y)
rootfinder = find_roots(spline)
optimafinder = find_optima(spline)
length(rootfinder.roots) == 1
length(optimafinder.optima) == 0
# But it has a point at which it has a value of four:
fourfinder = find_roots(spline; root_value = 4.0)
abs(evaluate(spline, fourfinder.roots[1]) - 4.0) < 1e-10
# and no points where it is negative four::
negfourfinder = find_roots(spline; root_value = -4.0)
length(negfourfinder.roots) == 0
fourfinder2 = find_roots(spline - 2.0; root_value = 2.0)
abs(fourfinder[:roots][1] - fourfinder2[:roots][1]) < eps()
# And if we strict domain to after 2.5 then we will not find any roots.
rootfinder22 = find_roots(spline; interval = (2.5,Inf))
length(rootfinder22.roots) == 0

# This has a root but no optima:
y = y .-2.0
spline2 = Schumaker(x,y)
rootfinder = find_roots(spline2)
optimafinder = find_optima(spline2)
length(rootfinder.roots) == 1
abs(rootfinder.roots[1] - 1.878) < 0.001
length(optimafinder.optima) == 0

y = (x .- 3).^2 .+ 6 # Should be an optima at x = 3. But no roots.
spline3 = Schumaker(x,y)
rootfinder = find_roots(spline3)
optimafinder = find_optima(spline3)
length(rootfinder.roots) == 0
length(optimafinder.optima) == 1
abs(optimafinder.optima[1] - 3.0) < 1e-2
optimafinder.optima_types[1] == :Minimum
