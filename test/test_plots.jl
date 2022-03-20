using Test
@testset "Test Plotting" begin
    using SchumakerSpline
    x = Array{Union{Missing,Float64}}(collect(0.0:0.01:2.0))
    y = Array{Union{Missing,Float64}}(sqrt.(x))
    y[5] = missing
    x[2] = missing
    s1 = Schumaker{Float64}(x, y)

    # This should plot something over [0,2]
    plt = plot(s1)
    # This should plot something over [0,1] but only two points because it is input as an array and this is the exact array of plotpoints.
    plt = plot(s1, [0.0,1.0])
    # Now it is a tuple and so there will be grid_len intermediate points put it (200 by default)
    plt = plot(s1, (0.0,1.0))
    # And we can plot the derivatives too
    plt = plot(s1, (0.0,1.0); derivs = true)
    # This should add a second spline just below with same derivaitves
    glt = plot(s1 -10.0, [0.5,1.0]; plot_options = (label = "shifted",), deriv_plot_options = (label = "shifted deriv1",), deriv2_plot_options = (label = "shifted deriv2",), plt = plt)
    # This should only plot the first spline.
    qlt = plot(s1 + 0.2, (0.0,1.0); plot_options = (label = "shifted",))


    ss = Array{Schumaker,1}([s1, s1+0.2, s1*0.8+0.2])
    plt  = plot(ss, [0.0,1.0])

    plt = plot(ss)
    plt2 = plot(0-2.0*s1 + 0.1, [0.0,1.0]; derivs = false, plot_options = (label = "new one",), plt = plt)
end
