module SchumakerSpline

    #using Plots
    using Dates

    include("SchumakerFunctions.jl")

    export Schumaker_ExtrapolationSchemes, Curve, Linear, Constant
    export Schumaker, evaluate, evaluate_integral
    include("roots_optima_intercepts.jl")
    export find_derivative_spline, find_roots, find_optima, get_crossover_in_interval, get_intersection_points
    include("map_to_shape.jl")
    export reshape_values
    include("splice_splines.jl")
    export splice_splines
    include("algebra.jl")
    export +,-,*,/
    include("plotting.jl")
    export plot
    include("higher_dimensions.jl")
    export Schumaker2d
end
