using SchumakerSpline
using Test

# Run tests


println("Test with Float64s")
include("Test_with_Floats.jl")
println("Test with Ints")
include("Test_with_Ints.jl")
println("Test with Dates")
include("Test_with_Dates.jl")
println("Test Intercepts")
include("test_intercepts.jl")
println("Test Splicing of Splines")
include("test_splice_splines.jl")
println("Test Plots")
include("test_plots.jl")
println("Test Extrapolation")
include("test_extrapolation.jl")
println("Test mapping to shape")
include("test_map_to_shape.jl")
println("Test Algebra")
include("test_algebra.jl")
println("Test 2d")
include("test_2d.jl")
