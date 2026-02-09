"""
test_if_intercept_in_interval(a1::Real,b1::Real,c1::Real,c2::Real,interval_width::Real)
    This tests if a spline could have passed over zero in a certain interval. The a1,b1,c1 are the coefficients of the spline. The two xs are for the left and right and c2 is the right hand level.
        Note that this function will not detect zeros that are precisely on the endpoints.
"""
function test_if_intercept_in_interval(a1::Real,b1::Real,c1::Real,c2::Real,interval_width::Real)
    if (sign(c1) == 0) || abs(sign(c1) - sign(c2)) > 1.5 return true end # If we cross the barrier then there is at least one intercept in interval.
    if sign(b1) == sign(2*a1*(interval_width)+b1) return false end # If we did not cross the barrier and the spline is monotonic then we did not cross
    # Now we have the case where the gradient switches sign within an interval but the sign of the endpoints did not change.
    # The easiest way to test will be to find the vertex of the parabola. See if it is within the interval and of a different sign to the endpoints.
    # We don't actually have to test if the vertex is in the interval however - it has to be for the gradient sign to have flipped.
    vertex_x = -b1/(2*a1) # Note that this is relative to x1.
    vertex_y = a1 * (vertex_x)^2 + b1*(vertex_x) + c1
    is_vertex_of_opposite_sign_in_y = abs(sign(c1) - sign(vertex_y)) > 0.5
    return is_vertex_of_opposite_sign_in_y
end

"""
    find_roots(spline::Schumaker{T}; root_value::Real = 0.0, interval::Tuple{<:Real,<:Real} = (spline.IntStarts_[1], spline.IntStarts_[length(spline.IntStarts_)])) where T<:Real
Finds roots - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find roots.
Here root_value can be set to get all points at which the function is equal to the root value. For instance if you want to find all points at which
the spline has a value of 1.0.
### Inputs
* `spline` - The spline you want to find the roots for.
* `root_value` - What level counts as a root.
* `interval` - What interval to explore for roots.
### Returns
* A `NamedTuple` describing all roots found together with the derivatives and second derivatives at that point.
"""
function find_roots(spline::Schumaker{T}; root_value::Real = 0.0, interval::Tuple{<:Real,<:Real} = (spline.IntStarts_[1], spline.IntStarts_[length(spline.IntStarts_)])) where T<:Real
    roots = T[]
    first_derivatives = T[]
    second_derivatives = T[]
    first_interval_start = searchsortedlast(spline.IntStarts_, interval[1])
    last_interval_start  = searchsortedlast(spline.IntStarts_, interval[2])
    len = length(spline.IntStarts_)
    go_from = max(1,first_interval_start)
    go_until = last_interval_start < len ? last_interval_start : len-1
    for i in go_from:go_until
        @inbounds a1 = spline.coefficient_matrix_[i,1]
        @inbounds b1 = spline.coefficient_matrix_[i,2]
        @inbounds c1 = spline.coefficient_matrix_[i,3] - root_value
        @inbounds c2 = spline.coefficient_matrix_[i+1,3] - root_value
        @inbounds interval_width = spline.IntStarts_[i+1] - spline.IntStarts_[i] + 1000*eps()
        if test_if_intercept_in_interval(a1,b1,c1,c2,interval_width)
            if abs(a1) > eps() # Is it quadratic
                det = sqrt(max(0.0, b1^2 - 4*a1*c1))
                r1 = (-b1 + det) / (2*a1)
                r2 = (-b1 - det) / (2*a1)
                left_root  = min(r1, r2)
                right_root = max(r1, r2)
                if (left_root >= 0) && (left_root <= interval_width)
                    @inbounds push!(roots, spline.IntStarts_[i] + left_root)
                    push!(first_derivatives, 2 * a1 * left_root + b1)
                    push!(second_derivatives, 2 * a1)
                end
                if (right_root >= 0) && (right_root <= interval_width)
                    @inbounds push!(roots, spline.IntStarts_[i] + right_root)
                    push!(first_derivatives, 2 * a1 * right_root + b1)
                    push!(second_derivatives, 2 * a1)
                end
            else # Is it linear?
                @inbounds new_root = spline.IntStarts_[i] - c1/b1
                if !((length(roots) > 0) && (abs(new_root - last(roots)) < 1e-5))
                    push!(roots, new_root)
                    push!(first_derivatives, b1)
                    push!(second_derivatives, zero(T))
                end
            end
        end
    end
    # Now adding on roots that occur after the end of the last interval.
    end_of_last_interval = spline.IntStarts_[len]
    if interval[2] >= end_of_last_interval
        @inbounds a = spline.coefficient_matrix_[len,1]
        @inbounds b = spline.coefficient_matrix_[len,2]
        @inbounds c = spline.coefficient_matrix_[len,3] - root_value
        if abs(a) > eps() # Is it quadratic
            root_determinant = sqrt(max(0.0, b^2 - 4*a*c))
            er1 = end_of_last_interval + (-b - root_determinant)/(2*a)
            er2 = end_of_last_interval + (-b + root_determinant)/(2*a)
            for er in (er1, er2)
                if er >= end_of_last_interval && er <= interval[2]
                    push!(roots, er)
                    push!(first_derivatives, 2 * a * (er - end_of_last_interval) + b)
                    push!(second_derivatives, 2 * a)
                end
            end
        elseif abs(b) > eps() # If it is linear.
            nr = -c/b + end_of_last_interval
            if nr >= end_of_last_interval && nr <= interval[2]
                push!(roots, nr)
                push!(first_derivatives, 2 * a * (nr - end_of_last_interval) + b)
                push!(second_derivatives, 2 * a)
            end
        end # We do nothing in the case that we have a constant - no chance of root.
    end
    # Sometimes if there are two roots within an interval and the endpoint of the interval is also here we get too many roots.
    # So here we get rid of stuff we don't want.
    if length(roots) == 0
        return (roots = roots, first_derivatives = first_derivatives, second_derivatives = second_derivatives)
    else
        nroots = length(roots)
        keep = trues(nroots)
        for i in 1:nroots
            if roots[i] < interval[1] || roots[i] > interval[2]
                keep[i] = false
            end
        end
        if nroots > 1
            for i in 1:(nroots-1)
                if abs(roots[i+1] - roots[i]) < 10000 * eps()
                    keep[i+1] = false
                end
            end
        end
        return (roots = roots[keep], first_derivatives = first_derivatives[keep], second_derivatives = second_derivatives[keep])
    end
end

"""
    find_optima(spline::Schumaker)
Finds optima - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find optima.
### Inputs
* `spline` - The spline you want to find optima for.
* `interval` - The interval over which you want to look for optima.
### Returns
* A NamedTuple containing the optima and the types of the optima (:Maximum or :Minimum)
"""
function find_optima(spline::Schumaker; interval::Tuple{<:Real,<:Real} = (spline.IntStarts_[1], spline.IntStarts_[length(spline.IntStarts_)]))
    deriv_spline = find_derivative_spline(spline)
    root_info = find_roots(deriv_spline; interval = interval)
    optima = root_info.roots
    optima_types =  Array{Symbol,1}(undef,length(optima))
    for i in 1:length(optima)
        if root_info.first_derivatives[i] > 1e-15
            optima_types[i] = :Minimum
        elseif root_info.first_derivatives[i] < -1e-15
            optima_types[i] = :Maximum
        else
            optima_types[i] = :SaddlePoint
        end
    end
    return (optima = optima, optima_types = optima_types)
end

## Finding intercepts
"""
    quadratic_formula_roots(a::Real,b::Real,c::Real)
A basic application of the textbook quadratic formula.
### Inputs
* `a` - The quadratic term
* `b` - The linear term
* `c` - The constant
### Returns
* A vector with the roots.
"""
function quadratic_formula_roots(a::Real,b::Real,c::Real)
    if abs(a) <= eps()
        if abs(b) <= eps()
            return Float64[]
        end
        return [-c/b]
    end
    discriminant = b^2 - 4*a*c
    determin = sqrt(max(0.0, discriminant))
    roots = [(-b + determin)/(2*a), (-b - determin)/(2*a)]
    return roots
end
"""
    get_crossover_in_interval(s1::Schumaker{T}, s2::Schumaker{R}, interval::Tuple{U,U}) where T<:Real where R<:Real where U<:Real
Finds the point at which two schumaker splines cross over each other within a single interval.
### Inputs
* `s1` - The first spline
* `s2` - The second spline
* `interval` - The interval you want to examine for crossovers.
### Returns
* A `Vector` describing crossover points.
"""
function get_crossover_in_interval(s1::Schumaker{T}, s2::Schumaker{R}, interval::Tuple{U,U}) where T<:Real where R<:Real where U<:Real
    # Getting the coefficients for the first spline.
    i = searchsortedlast(s1.IntStarts_, interval[1])
    @inbounds start1 = s1.IntStarts_[i]
    @inbounds a1 = s1.coefficient_matrix_[i,1]
    @inbounds b1 = s1.coefficient_matrix_[i,2]
    @inbounds c1 = s1.coefficient_matrix_[i,3]
    # Getting the coefficients for the second spline.
    j = searchsortedlast(s2.IntStarts_, interval[1])
    @inbounds start2 = s2.IntStarts_[j]
    @inbounds a2 = s2.coefficient_matrix_[j,1]
    @inbounds b2 = s2.coefficient_matrix_[j,2]
    @inbounds c2 = s2.coefficient_matrix_[j,3]
    # Get implied coefficients for the s1 - s2 quadratic. Pretty simple algebra gets this.
    # As a helper we define G = start2 - start1. We define A,B,C as coefficients of s1-s2.
    # The final spline is in terms of (x-start1).
    G = start1 - start2
    A = a1 - a2
    B = b1 - b2 - 2*a2*G
    C = c1 - c2 - a2*(G^2) - b2*G
    # Now we need to use quadratic formula to get the roots and pick the root in the interval.
    roots = quadratic_formula_roots(A,B,C) .+ start1
    roots_in_interval = roots[(roots .>= interval[1]-10*eps()) .& (roots .<= interval[2]+10*eps())]
    return roots_in_interval
end

"""
    get_intersection_points(s1::Schumaker{T}, s2::Schumaker{R}) where T<:Real where R<:Real
This funds the coordinates of the point at which spline s1 intercepts spline s2.
### Inputs
* `s1` - The first spline
* `s2` - The second spline
### Returns
* Locations of any crossover points.
"""
function get_intersection_points(s1::Schumaker{T}, s2::Schumaker{R}) where T<:Real where R<:Real
    # What x locations to loop over
    all_starts = sort(unique(vcat(s1.IntStarts_, s2.IntStarts_)))
    start_of_overlap = max(s1.IntStarts_[1], s2.IntStarts_[1])
    overlap_starts = all_starts[all_starts .>= start_of_overlap]
    # Getting a container to return results
    promo_type = promote_type(T,R)
    locations_of_crossovers = Array{promo_type,1}()
    # For the first part what function is higher.
    last_one_greater = evaluate(s1, overlap_starts[1]) > evaluate(s2, overlap_starts[1])
    for i in 2:length(overlap_starts)
        start = overlap_starts[i]
        val_1 = evaluate(s1, start)
        val_2 = evaluate(s2, start)
        # Need to take into account the ordering and record when it flips
        one_greater = val_1 > val_2
        if one_greater != last_one_greater
            interval = (overlap_starts[i-1], overlap_starts[i])
            crossover = get_crossover_in_interval(s1, s2, interval)
            if length(crossover) != 1
                error("Only one crossover expected in interval from a continuous spline.")
            end
            push!(locations_of_crossovers, crossover[1])
        end
        last_one_greater = one_greater
    end
    return locations_of_crossovers
end
