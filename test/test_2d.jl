using Test
@testset "Test two dimensional" begin
    using SchumakerSpline

    gridx = collect(1:10)
    gridy = collect(1:10)
    grid = gridx * gridy'

    schum = Schumaker2d(gridx, gridy, grid; extrapolation = (Linear, Linear))

    @test abs(schum(5,5) - 25) < 100 * eps()
    @test abs(schum(5.5,5.5) - 30.25) < 100 * eps()
end
