module SchumakerSpline

using Dates

include("SchumakerFunctions.jl")

export Schumaker_ExtrapolationSchemes, Curve, Linear, Constant
export Schumaker, evaluate, evaluate_integral

end
