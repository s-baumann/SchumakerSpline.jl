using SchumakerSpline
using Base.Test

# Run tests

tic()
println("Test with Float64s")
@time @test include("Test_with_Float64s.jl")
toc()
tic()
println("Test with Ints")
@time @test include("Test_with_Ints.jl")
toc()
tic()
println("Test with Dates")
@time @test include("Test_with_Dates.jl")
toc()
