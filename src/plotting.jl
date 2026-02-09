const _CM_TO_PX = 37.795275591  # 96 DPI

"""
    plot(s1::Schumaker, interval::Tuple{R,R}; derivs=false, grid_len=200, label="Spline",
         deriv_label="1st derivative", deriv2_label="2nd derivative",
         width_in_cm=20, height_in_cm=12) where R<:Real
    plot(s1::Schumaker, grid::AbstractArray{R,1}; derivs=false, label="Spline",
         deriv_label="1st derivative", deriv2_label="2nd derivative",
         width_in_cm=20, height_in_cm=12) where R<:Real

### Inputs
* `s1` - The Schumaker spline to chart
* `interval` - The interval over which to chart it.
* `grid_len` - The number of grid points to be used in plotting.
* `grid` - The grid. If used this is instead of the `interval` and `grid_len`
* `derivs` - Should the derivative splines also be plotted
* `label` - Label for the spline curve
* `deriv_label` - Label for the first derivative curve
* `deriv2_label` - Label for the second derivative curve
* `width_in_cm` - Width of the plot in centimeters (default 20)
* `height_in_cm` - Height of the plot in centimeters (default 12)
### Returns
* A VegaLite plot specification
"""
function plot(s1::Schumaker, interval::Tuple{R,R} = (s1.IntStarts_[1], s1.IntStarts_[length(s1.IntStarts_)]); derivs::Bool = false, grid_len::Integer = 200,
              label::String = "Spline", deriv_label::String = "1st derivative", deriv2_label::String = "2nd derivative",
              width_in_cm::Real = 20, height_in_cm::Real = 12) where R<:Real
    grid = collect(range(interval[1], interval[2], length=grid_len))
    return plot(s1, grid; derivs = derivs, label = label, deriv_label = deriv_label, deriv2_label = deriv2_label,
                width_in_cm = width_in_cm, height_in_cm = height_in_cm)
end
function plot(s1::Schumaker, grid::AbstractArray{R,1}; derivs::Bool = false,
              label::String = "Spline", deriv_label::String = "1st derivative", deriv2_label::String = "2nd derivative",
              width_in_cm::Real = 20, height_in_cm::Real = 12) where R<:Real
    evals = s1.(grid)
    data = [(x = grid[i], y = evals[i], series = label) for i in eachindex(grid)]
    if derivs
        evals_1 = evaluate.(s1, grid, 1)
        evals_2 = evaluate.(s1, grid, 2)
        append!(data, [(x = grid[i], y = evals_1[i], series = deriv_label) for i in eachindex(grid)])
        append!(data, [(x = grid[i], y = evals_2[i], series = deriv2_label) for i in eachindex(grid)])
    end
    w = round(Int, width_in_cm * _CM_TO_PX)
    h = round(Int, height_in_cm * _CM_TO_PX)
    return data |> @vlplot(:line, x = :x, y = :y, color = :series, width = w, height = h)
end


"""
    plot(ss::Vector{<:Schumaker}, interval::Tuple{R,R}; derivs=false, grid_len=200,
         labels=missing, width_in_cm=20, height_in_cm=12) where R<:Real
    plot(ss::Vector{<:Schumaker}, grid::AbstractArray{R,1}; derivs=false,
         labels=missing, width_in_cm=20, height_in_cm=12) where R<:Real
### Inputs
* `ss` - a vector of Schumaker splines to chart
* `interval` - The interval over which to chart it.
* `grid_len` - The number of grid points to be used in plotting.
* `grid` - The grid. If used this is instead of the `interval` and `grid_len`
* `derivs` - Should the derivative splines also be plotted
* `labels` - Labels for the spline curves. Defaults to "Spline 1", "Spline 2", etc.
* `width_in_cm` - Width of the plot in centimeters (default 20)
* `height_in_cm` - Height of the plot in centimeters (default 12)
### Returns
* A VegaLite plot specification
"""
function plot(ss::Vector{<:Schumaker}, interval::Tuple{R,R} = (ss[1].IntStarts_[1], ss[1].IntStarts_[length(ss[1].IntStarts_)]); derivs::Bool = false, grid_len::Integer = 200,
              labels::Union{Vector{String},Missing} = missing, width_in_cm::Real = 20, height_in_cm::Real = 12) where R<:Real
    grid = collect(range(interval[1], interval[2], length=grid_len))
    return plot(ss, grid; derivs = derivs, labels = labels, width_in_cm = width_in_cm, height_in_cm = height_in_cm)
end

function plot(ss::Vector{<:Schumaker}, grid::AbstractArray{R,1}; derivs::Bool = false,
              labels::Union{Vector{String},Missing} = missing, width_in_cm::Real = 20, height_in_cm::Real = 12) where R<:Real
    data = NamedTuple{(:x, :y, :series), Tuple{Float64, Float64, String}}[]
    for i in eachindex(ss)
        lbl = ismissing(labels) ? "Spline $i" : labels[i]
        evals = ss[i].(grid)
        append!(data, [(x = grid[j], y = evals[j], series = lbl) for j in eachindex(grid)])
        if derivs
            evals_1 = evaluate.(ss[i], grid, 1)
            evals_2 = evaluate.(ss[i], grid, 2)
            append!(data, [(x = grid[j], y = evals_1[j], series = "$lbl - 1st deriv") for j in eachindex(grid)])
            append!(data, [(x = grid[j], y = evals_2[j], series = "$lbl - 2nd deriv") for j in eachindex(grid)])
        end
    end
    w = round(Int, width_in_cm * _CM_TO_PX)
    h = round(Int, height_in_cm * _CM_TO_PX)
    return data |> @vlplot(:line, x = :x, y = :y, color = :series, width = w, height = h)
end
