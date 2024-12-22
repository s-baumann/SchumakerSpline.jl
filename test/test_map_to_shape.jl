using Test
@testset "Test mapping to a shape" begin
    using SchumakerSpline

    xvals = [1,2,3,4,5,6]
    shape_map(x,upper,lower) = min(max(x,lower), upper)

    # Increasing
    yvals = [1.0, 2.0, 2.9, 3.9, 2.3, 1.4]
    vals2 = reshape_values(xvals, yvals; increasing = true, concave = true, shape_map = shape_map)
    @test all(abs.([1.0, 2.0, 2.9, 3.8, 3.8, 3.8] .- vals2) .< 10*eps()) # Changes

    # Increasing and concave
    yvals = [1.0, 2.0, 2.9, 3.9, 4.3, 4.4]
    vals2 = reshape_values(xvals, yvals; increasing = true, concave = true, shape_map = shape_map)
    @test all(abs.([1.0, 2.0, 2.9, 3.8, 4.3, 4.4] .- vals2) .< 10*eps()) # Changes
    yvals = sqrt.(xvals)
    vals2 = reshape_values(xvals, yvals; increasing = true, concave = true, shape_map = shape_map)
    @test all(abs.(vals2 .- yvals) .< 10*eps()) # No changes as input already increasing concave

    # Increasing and convex
    yvals = [1.0, 2.0, 3.1, 3.9, 4.3, 4.4]
    vals2 = reshape_values(xvals, yvals; increasing = true, concave = false, shape_map = shape_map)
    @test all(abs.([1.0, 2.0, 3.1, 4.2, 5.3, 6.4] .- vals2) .< 10*eps()) #changes
    yvals = (xvals) .^ 2
    vals2 = reshape_values(xvals, yvals; increasing = true, concave = false, shape_map = shape_map)
    @test all(abs.(vals2 .- yvals) .< 10*eps()) # No changes as input already increasing concave


    # Decreasing
    yvals = -[1.0, 2.0, 2.9, 3.9, 2.3, 1.4]
    vals2 = reshape_values(xvals, yvals; increasing = false, concave = false, shape_map = shape_map)
    @test all(abs.([-1.0, -2.0, -2.9, -3.8, -3.8, -3.8] .- vals2) .< 10*eps()) # Changes

    # Decreasing and convex
    yvals = -[1.0, 2.0, 2.9, 3.9, 4.3, 4.4]
    vals2 = reshape_values(xvals, yvals; increasing = false, concave = false, shape_map = shape_map)
    @test all(abs.(-[1.0, 2.0, 2.9, 3.8, 4.3, 4.4] .- vals2) .< 10*eps()) # Changes
    yvals = -sqrt.(xvals)
    vals2 = reshape_values(xvals, yvals; increasing = false, concave = false, shape_map = shape_map)
    @test all(abs.(vals2 .- yvals) .< 10*eps()) # No changes as input already increasing concave

    # Decreasing and concave
    yvals = -[1.0, 2.0, 3.1, 3.9, 4.3, 4.4]
    vals2 = reshape_values(xvals, yvals; increasing = false, concave = true, shape_map = shape_map)
    @test all(abs.(-[1.0, 2.0, 3.1, 4.2, 5.3, 6.4] .- vals2) .< 10*eps()) #changes
    yvals = -(xvals) .^ 2
    vals2 = reshape_values(xvals, yvals; increasing = false, concave = true, shape_map = shape_map)
    @test all(abs.(vals2 .- yvals) .< 10*eps()) # No changes as input already increasing concave
end
