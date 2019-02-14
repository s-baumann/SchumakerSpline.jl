module SchumakerSpline

using Dates

include("SchumakerFunctions.jl")

export Schumaker_ExtrapolationSchemes, Curve, Linear, Constant
export Schumaker, evaluate, evaluate_integral
export find_derivative_spline, find_roots, find_optima
end
