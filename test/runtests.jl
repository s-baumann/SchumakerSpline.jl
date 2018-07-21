using SchumakerSpline
using Base.Test

# Run tests

tic()
println("Test with Float64s")
@time @test include("Test_with_Float64s.jl")
toc()
tic()
println("Test with Int64s")
@time @test include("Test_with_Int64s.jl")
toc()
tic()
println("Test with Dates")
@time @test include("Test_with_Dates.jl")
toc()
