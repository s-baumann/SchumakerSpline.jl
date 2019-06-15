"""
find_root(spline::Schumaker; root_value::T = 0.0)
    Finds roots - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find roots.
    Here root_value can be set to get all points at which the function is equal to the root value. For instance if you want to find all points at which
    the spline has a value of 1.0.
"""
function find_roots(spline::Schumaker; root_value::T = 0.0) where T<:Real
    roots = Array{Float64,1}(undef,0)
    first_derivatives = Array{Float64,1}(undef,0)
    second_derivatives = Array{Float64,1}(undef,0)
    len = length(spline.IntStarts_)
    constants = spline.coefficient_matrix_[:,3]
    if len < 2
        return (roots = roots, first_derivatives = first_derivatives, second_derivatives = second_derivatives)
    else
        for i in 1:(len-1)
            if abs(sign(constants[i] - root_value) - sign(constants[i+1] - root_value)) > 0.5
                a = spline.coefficient_matrix_[i,1]
                b = spline.coefficient_matrix_[i,2]
                c = spline.coefficient_matrix_[i,3] - root_value
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
    end
    return (roots = roots, first_derivatives = first_derivatives, second_derivatives = second_derivatives)
end

"""
find_optima(spline::Schumaker)
Finds optima - This is handy because in many applications schumaker splines are monotonic and globally concave/convex and so it is easy to find optima.

"""
function find_optima(spline::Schumaker)
    deriv_spline = find_derivative_spline(spline)
    root_info = find_roots(deriv_spline)
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
