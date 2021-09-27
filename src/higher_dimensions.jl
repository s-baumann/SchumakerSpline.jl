
"""
Experimental: This uses a combination of Schumaker Splines to cover a 2 dimensional space.
We have a grid of splines. We first evaluate in one dimension (which leads us to a point between two adjacent splines).
Then we evaluate each of the two splines and interpolate.
### Members
* `IntStarts_` - A vector with the coordinates of each schumaker spline.
* `schumakers` - A vector of schumaker splines.
"""
struct Schumaker2d{T<:AbstractFloat}
    IntStarts_::Array{T,1}
    schumakers::Array{Schumaker{T},1}
    function Schumaker2d(x_row::Vector{R}, x_col::Vector{Y}, ygrid::Array{T,2}; bycol=true,
             extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve)) where T<:Real where R<:Real where Y<:Real
        pro = promote_type(promote_type(T,Y),R)
        shape = size(ygrid)
        dimension = 1 + Integer(bycol)
        schums = Array{Schumaker{pro},1}(undef, shape[dimension])
        for i in 1:shape[dimension]
           ys = bycol ?  ygrid[:,i] : ygrid[i,:]
           xs = bycol ? x_row : x_col
           schums[i] = Schumaker(xs, ys; extrapolation = extrapolation)
        end
        return new{pro}(bycol ? x_col : x_row,schums)
    end
    function Schumaker2d(x_row::Vector{R}, x_col::Vector{Y}, ygrid::Array{T,2}; bycol=true,
             extrapolation::Tuple{Schumaker_ExtrapolationSchemes,Schumaker_ExtrapolationSchemes} = (Curve,Curve)) where T<:Integer where R<:Integer where Y<:Integer
        fl = typeof(AbstractFloat(1.0))
        return Schumaker2d( fl.(x_row), fl.(x_col), fl.(ygrid); bycol = bycol, extrapolation = extrapolation  )
    end
end

"""
    evaluate(spline::Schumaker2d, p1::Real, p2::Real)
### Inputs
* `spline` - A Schumaker2d
* `p1` - The coordinate in the first dimension.
* `p2` - The coordinate in the second dimension.
### Returns
* A scalar
"""
function evaluate(spline::Schumaker2d, p1::Real, p2::Real)
    distances   = abs.(spline.IntStarts_ .- p1)
    # Direct hits
    direct_hits = findall(distances .< 100*eps())
    if length(direct_hits) > 0
        return spline.schumakers[direct_hits[1]](p2)
    end
    closest = findall(distances .< sort(distances)[2] + 100*eps())[[1,2]]
    dists = distances[closest]
    ys = map(s -> s(p2), spline.schumakers[closest])
    y = sum(ys .* ( reverse(dists)  ) ./ sum(dists))
    return y
end
function (s::Schumaker2d)(p1::Real, p2::Real)
    return evaluate(s, p1, p2)
end
