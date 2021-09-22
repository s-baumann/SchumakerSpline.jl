
shape_map(x,upper,lower) = min(max(x,lower), upper)


"""
    reshape_values(xvals::Vector{<:Real}, yvals::Vector{<:Real}; increasing::Bool = true,
                       concave::Bool = true, shape_map::Function = shape_map)
This reshapes a vector of yvalues. For instance if we are doing fixed point acceleration that
should result in a monotonic concave function (ie the consumption smoothing problem from the documentation examples)
then we may end up with occasional non monotonic/concave values due to a dodgy optimiser or some other numerical issue.
So we can use reshape_values to adjust the values that cannot be true.
### Inputs
* `xvals` - A vector of x coordinates.
* `yvals` - A vector of y coordinates
* `increasing` - Should the y values be increasing. If false then they must be decreasing
* `concave` - Should the y values be concave. If false then they must be convex
* `shape_map` - A function used to adjust values to be increasing-concave (or whatever settings)
### Returns
* An updated vector of y values.
"""
function reshape_values(xvals::Vector{<:Real}, yvals::Vector{<:Real}; increasing::Bool = true,
                        concave::Bool = true, shape_map::Function = shape_map)
    lenlen = length(yvals)
    new_shape_vec = Array{Float64,1}(undef, lenlen)
    if lenlen > 0 new_shape_vec[1] = yvals[1] end
    if lenlen > 1
        upper = increasing ? Inf : new_shape_vec[1]
        lower = increasing ? new_shape_vec[1] : -Inf
        new_shape_vec[2] = shape_map(yvals[2],upper, lower)
    end
    for i in 3:lenlen
        stepp = xvals[i] - xvals[i-1]
        previous_grad = (new_shape_vec[i-1]-new_shape_vec[i-2])/(xvals[i-1] - xvals[i-2])
        upper =  new_shape_vec[i-1] + (increasing ? (concave ? previous_grad * stepp : Inf) : (concave ? previous_grad * stepp : 0))
        lower =  new_shape_vec[i-1] + (increasing ? (concave ? 0 :  previous_grad * stepp) : (concave ? -Inf : previous_grad * stepp))
        new_shape_vec[i] = shape_map(yvals[i],upper, lower)
    end
    return new_shape_vec
end
