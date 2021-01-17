
function reshape_values(xvals::Vector{<:Real}, yvals::Vector{<:Real}; increasing::Bool = true,
                        concave::Bool = true, shape_map::Function = min(max(x,lower), upper))
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
