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

# This is a historical bug - It did not find the root because it happened in the interval right before 4.0
spline = Schumaker{Float64}([0.0, 0.166667, 0.25, 0.277597, 0.5, 0.581563, 0.75, 0.974241, 1.0, 95.1966, 100.0, 116.344, 200.0, 233.333], [1.92913 -9.4624 9.91755; 7.71653 -23.3254 8.39407; 40.1234 -101.466 6.50388; 0.617789 -3.15302 3.73426; 1.84611 -5.68404 3.06358; 0.432882 -2.23154 2.61226; 0.3396 -1.88818 2.24866;
                               25.7353 -51.7173 1.84233; 2.00801e-5 -0.183652 0.527206; 0.00772226 0.288784 -16.594; 0.00225409 0.0127389 -15.0287; 8.60438e-5 -0.0760628 -14.2184; 0.000216017 -0.0446771 -19.9793; 5.40043e-5 -0.0550539 -21.2285])
interval = (1e-14, 4.0)
root_value = 0.0
optima = find_roots(spline; root_value = root_value, interval = interval)
length(optima.roots) == 1
abs(optima.roots[1] - 3.875) < 0.005


spline = Schumaker{Float64}([0.0, 0.166667, 0.25, 0.277597, 0.5, 0.581563, 0.75, 0.913543, 1.0, 1.65249, 2.0, 2.33612, 3.0, 3.5334, 4.0, 4.66667], [-0.645627 2.21521 0.0; -2.58251 2.0 0.351267; -13.4282 1.56958 0.5; -0.206757 0.828427 0.533089; -0.617841 0.736461 0.707107; -0.144873 0.635674 0.763065;
                             -0.155837 0.58687 0.866025; -0.557609 0.535898 0.957836; -0.0193613 0.43948 1.0; -0.068257 0.414214 1.27851; -0.072796 0.366774 1.41421; -0.0186602 0.317837 1.52927; -0.0235395 0.293061 1.73205; -0.030761 0.267949 1.88167; -0.0108734 0.239243 2.0; -0.00271836 0.224745 2.15466])
root_value = 2.0
optima = find_roots(spline; root_value = root_value)
length(optima.roots) == 1
abs(optima.roots[1] - 4.0) < 1e-10
