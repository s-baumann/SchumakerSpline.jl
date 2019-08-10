"""
find_root(spline::Schumaker; root_value::T = 0.0)
    Finds roots - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find roots.
    Here root_value can be set to get all points at which the function is equal to the root value. For instance if you want to find all points at which
    the spline has a value of 1.0.
"""
function find_roots(spline::Schumaker{T}; root_value::Real = 0.0, interval::Tuple{<:Real,<:Real} = (spline.IntStarts_[1], spline.IntStarts_[length(spline.IntStarts_)])) where T<:Real
    roots = Array{T,1}(undef,0)
    first_derivatives = Array{T,1}(undef,0)
    second_derivatives = Array{T,1}(undef,0)
    first_interval_start = searchsortedlast(spline.IntStarts_, interval[1])
    last_interval_start  = searchsortedlast(spline.IntStarts_, interval[2])
    len = length(spline.IntStarts_)
    go_until = last_interval_start < len ? last_interval_start : len-1
    constants = spline.coefficient_matrix_[:,3]
    constants_minus_root = constants .- root_value
    for i in first_interval_start:go_until
        if abs(constants_minus_root[i]) < eps()
            a = spline.coefficient_matrix_[i,1]
            b = spline.coefficient_matrix_[i,2]
            root = spline.IntStarts_[i]
            append!(roots, root)
            append!(first_derivatives, 2 * a * root + b)
            append!(second_derivatives, 2 * a)
            continue # We don't need to do the next bit.
        end
        if abs(sign(constants_minus_root[i]) - sign(constants_minus_root[i+1])) > 1.5 # 1.5 because if we have exact equalities we dont want to consider them except in the above code block.
            a = spline.coefficient_matrix_[i,1]
            b = spline.coefficient_matrix_[i,2]
            c = constants_minus_root[i]
            if abs(a) > 1e-13 # Is it quadratic
                det = sqrt(b^2 - 4*a*c)
                left_root  = (-b + det) / (2*a) # The x coordinate here is relative to spline.IntStarts_[i]. We want the smallest one that is to the right (ie positive)
                right_root = (-b - det) / (2*a)
                if (left_root > 1e-15) && (left_root < spline.IntStarts_[i+1] - spline.IntStarts_[i] +  1e-15)
                    append!(roots, left_root + spline.IntStarts_[i])
                    append!(first_derivatives, 2 * a * left_root + b)
                    append!(second_derivatives, 2 * a)
                elseif (right_root > 1e-15) && (right_root < spline.IntStarts_[i+1] - spline.IntStarts_[i] +  1e-15)
                    append!(roots, right_root + spline.IntStarts_[i])
                    append!(first_derivatives, 2 * a * right_root + b)
                    append!(second_derivatives, 2 * a)
                end
            else # Is it linear? Note it cannot be constant or else it could not have jumped past zero in the interval.
                new_root = spline.IntStarts_[i] - c/b
                if !((length(roots) > 0) && (abs(new_root - last(roots)) < 1e-5))
                    append!(roots, spline.IntStarts_[i] - c/b)
                    append!(first_derivatives, b)
                    append!(second_derivatives, 0.0)
                end
            end
        end
    end
    return (roots = roots, first_derivatives = first_derivatives, second_derivatives = second_derivatives)
end

"""
find_optima(spline::Schumaker)
Finds optima - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find optima.

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
            optima_types = :Maximum
        else
            optima_types = :SaddlePoint
        end
    end
    return (optima = optima, optima_types = optima_types)
end

## Finding intercepts
"""
    quadratic_formula_roots(a,b,c)
A basic application of the textbook quadratic formula.
"""
function quadratic_formula_roots(a,b,c)
    determin = sqrt(b^2 - 4*a*c)
    roots = [(-b + determin)/(2*a), (-b - determin)/(2*a)]
    return(roots)
end
"""
    get_crossover_in_interval(s1::Schumaker{T}, s2::Schumaker{R}, interval::Tuple{U,U}) where T<:Real where R<:Real where U<:Real
Finds the point at which two schumaker splines cross over each other within a single interval. This is not exported.

"""

function get_crossover_in_interval(s1::Schumaker{T}, s2::Schumaker{R}, interval::Tuple{U,U}) where T<:Real where R<:Real where U<:Real
    # Getting the coefficients for the first spline.
    i = searchsortedlast(s1.IntStarts_, interval[1])
    start1 = s1.IntStarts_[i]
    a1,b1,c1 = Tuple(s1.coefficient_matrix_[i,:])
    # Getting the coefficients for the second spline.
    j = searchsortedlast(s2.IntStarts_, interval[1])
    start2 = s2.IntStarts_[j]
    a2,b2,c2 = Tuple(s2.coefficient_matrix_[j,:])
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
    return(roots_in_interval)
end


function get_intersection_points(s1::Schumaker{T}, s2::Schumaker{R}) where T<:Real where R<:Real
    # What x locations to loop over
    all_starts = sort(unique(vcat(s1.IntStarts_, s2.IntStarts_)))
    start_of_overlap = maximum([minimum(s1.IntStarts_), minimum(s2.IntStarts_)])
    overlap_starts = all_starts[all_starts .> start_of_overlap]
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
            interval = Tuple([overlap_starts[i-1], overlap_starts[i]])
            crossover = get_crossover_in_interval(s1, s2, interval)
            @assert length(crossover) == 1 "Only one crossover expected in interval from a continuous spline."
            push!(locations_of_crossovers, crossover[1])
        end
        last_one_greater = one_greater
    end
    return locations_of_crossovers
end
