using SchumakerSpline
using Base.Test

# Run tests

tic()
println("Second Derivative Test")
@time @test include("SecondDerivativeTest.jl")
toc()
