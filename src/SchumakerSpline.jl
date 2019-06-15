module SchumakerSpline

using Dates

include("SchumakerFunctions.jl")

export Schumaker_ExtrapolationSchemes, Curve, Linear, Constant
export Schumaker, evaluate, evaluate_integral
include("roots_optima_intercepts.jl")
export find_derivative_spline, find_roots, find_optima, get_intersection_points
include("splice_splines.jl")
export splice_splines
end
