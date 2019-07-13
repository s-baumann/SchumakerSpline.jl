import Plots.plot

function plot(s1::Schumaker, interval::AbstractArray{R,1} = [s1.IntStarts_[1], 2*s1.IntStarts_[length(s1.IntStarts_)] - s1.IntStarts_[length(s1.IntStarts_)-1]]; # This default plots over the IntStarts.
              derivs::Bool = false, grid_len::Integer = 20, plot_options::NamedTuple = (label = "Spline",), deriv_plot_options::NamedTuple = (label = "Spline - 1st deriv",), deriv2_plot_options::NamedTuple = (label = "Spline - 2nd deriv",), plt = missing) where R<:Real
    if length(interval) != 2 error("The interval parameter must have two entries for the start and endpoint of the segment.") end
        interval_as_tuple = Tuple{R,R}([interval[1], interval[2]])
    return plot(s1,interval_as_tuple; derivs = derivs, grid_len = grid_len, plot_options = plot_options, deriv_plot_options = deriv_plot_options, deriv2_plot_options = deriv2_plot_options, plt = plt)
end
function plot(s1::Schumaker, interval::Tuple{R,R}; derivs::Bool = false, grid_len::Integer = 20, plot_options::NamedTuple = (label = "Spline",), deriv_plot_options::NamedTuple = (label = "Spline - 1st deriv",),
              deriv2_plot_options::NamedTuple = (label = "Spline - 2nd deriv",), plt = missing) where R<:Real
    grid = collect(range(interval[1], interval[2], length=grid_len))
    evals = evaluate.(s1,grid, 0)
    plt = ismissing(plt) ? plot() : plt
    plt = plot!(plt, grid, evals; plot_options...)
    if derivs
        evals_1 = evaluate.(s1,grid, 1)
        plt = plot!(plt, grid, evals_1; deriv_plot_options...)
        evals_2 = evaluate.(s1,grid, 2)
        plt = plot!(plt, grid, evals_2; deriv2_plot_options...)
        return plt
    else
        return plt
    end
end

function plot(ss::AbstractArray{Schumaker,1}, interval::AbstractArray{R,1} = [ss[1].IntStarts_[1], 2*ss[1].IntStarts_[length(ss[1].IntStarts_)] - ss[1].IntStarts_[length(ss[1].IntStarts_)-1]]; # This default plots over the IntStarts of the first spline.
    derivs::Bool = false, grid_len::Integer = 20, plot_options::Union{AbstractArray{Tuple,1},Missing} = missing, deriv_plot_options::Union{AbstractArray{Tuple,1},Missing} = missing,
    deriv2_plot_options::Union{AbstractArray{Tuple,1},Missing} = missing, plt = missing) where R<:Real
    if length(interval) != 2 error("The interval parameter must have two entries for the start and endpoint of the segment.") end
    interval_as_tuple = Tuple{R,R}([interval[1], interval[2]])
    return plot(ss,interval_as_tuple; derivs = derivs, grid_len = grid_len, plot_options = plot_options, deriv_plot_options = deriv_plot_options, deriv2_plot_options = deriv2_plot_options, plt = plt)
end

function plot(ss::AbstractArray{Schumaker,1}, interval::Tuple{R,R}; derivs::Bool = false, grid_len::Integer = 20, plot_options::Union{AbstractArray{Tuple,1},Missing} = missing, deriv_plot_options::Union{AbstractArray{Tuple,1},Missing} = missing,
              deriv2_plot_options::Union{AbstractArray{Tuple,1},Missing} = missing, plt = missing)  where R<:Real
    plt = ismissing(plt) ? plot() : plt
    for i in 1:length(ss)
        plot_options_i = ismissing(plot_options) ? (label = string("Spline ", i),) : plot_options[i]
        deriv_plot_options_i = ismissing(deriv_plot_options) ? (label = string("Spline ", i, " - 1st deriv"),) : deriv_plot_options[i]
        deriv2_plot_options_i = ismissing(deriv2_plot_options) ? (label = string("Spline ", i, " - 2nd deriv"),) : deriv2_plot_options[i]
        plt = plot(ss[i], interval; derivs = derivs, grid_len = grid_len, plot_options = plot_options_i, deriv_plot_options = deriv_plot_options_i, deriv2_plot_options = deriv2_plot_options_i, plt = plt)
    end
    return plt
end
