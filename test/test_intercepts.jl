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
# And if we strict domain to after 2.5 then we will find one root at 103ish
rootfinder22 = find_roots(spline; interval = (2.5,Inf))
spline(rootfinder22.roots[1]) < 1e-10

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

# This covers the case where there is an intercept but it comes after the last value of IntStarts but within the interval.
spline = Schumaker{Float64}([0.0, 0.166667, 0.25, 0.277597, 0.5, 0.581563, 0.75, 0.913543, 1.0, 1.75963, 2.0, 2.66667], [1.92913 -10.7562 11.8412; 7.71653 -28.5007 10.1021; 40.1234 -128.376 7.78064;
                          0.617789 -3.56736 4.2684; 1.84611 -6.92217 3.50557; 0.432882 -2.52187 2.95326; 0.465641 -2.46738 2.54076; 1.66614 -5.5315 2.1497; 0.0496917 -1.00572 1.68391; 0.496301 -1.87488 0.948612;
                          0.0929378 -0.847743 0.526628; 0.0232345 -0.618539 0.00277145])
interval = (1e-14, 4.0)
optima = find_roots(spline; interval = interval)
length(optima.roots) == 1
spline(optima.roots[1]) < 1e-10
# In this case the root is in the last interval which is linear rather than quadratic.
spline = Schumaker{Float64}([-1.0e-10, 0.0, 0.02, 0.03, 0.0460546, 0.1, 0.170253, 0.3, 0.361175, 0.5, 0.598784, 0.75, 0.881918, 1.0, 1.27089, 1.5, 1.8286, 2.0, 2.66667, 4.0], [-0.0 0.0 1.74211; 66.2224 -229.519 34.2346; 264.889 -872.979 29.6707; 182.141 -599.117 20.9674; 16.1324 -56.256 11.3957; 6.4478 -23.5354 8.40789; 1.89036 -8.26419 6.78629; 4.79607 -16.6813 5.74585; 0.931297 -4.40443 4.74333; 1.12989 -4.72817 4.14983; 0.482179 -2.64836 3.69379; 0.577267 -2.76551 3.30434; 0.720484 -2.99668 2.94957; 0.162497 -1.37639 2.60576; 0.227166 -1.44474 2.24484; 0.162381 -1.19387 1.92576; 0.59688 -1.97635 1.55098; 0.0393714 -0.694514 1.22978; 0.00984285 -0.598087 0.784266; -0.0 -0.583443 0.0121747])
interval = (1e-14, 5.0)
optima = find_roots(spline; interval = interval)
length(optima.roots) == 1
spline(optima.roots[1]) < 1e-10
