using Test
@testset "Regression tests" begin
    using SchumakerSpline

    # These tests pin exact numerical output to detect unintended changes
    # from refactoring or efficiency improvements. All expected values were
    # generated from the package at commit 1991efa (Julia 1.12).

    # ------------------------------------------------------------------
    # 1. Float spline: log(x) + sqrt(x), 1000 points on [1,6]
    # ------------------------------------------------------------------
    @testset "Float spline evaluate & integral" begin
        x = collect(range(1.0, stop=6.0, length=1000))
        y = log.(x) .+ sqrt.(x)
        sp = Schumaker(x, y)

        # Evaluation at on-grid and off-grid points (value, 1st & 2nd derivative)
        @test evaluate(sp, 1.0, 0)  == 1.0
        @test evaluate(sp, 1.0, 1)  == 1.4984288485329673
        @test evaluate(sp, 1.0, 2)  == -0.46366500258369364
        @test evaluate(sp, 1.5, 0)  == 1.6302099777999535
        @test evaluate(sp, 1.5, 1)  == 1.0749176285913955
        @test evaluate(sp, 1.5, 2)  == -0.5775096109828065
        @test evaluate(sp, 2.3, 0)  == 2.3494842109951404
        @test evaluate(sp, 2.3, 1)  == 0.7644728540519802
        @test evaluate(sp, 2.3, 2)  == -0.25974335141659244
        @test evaluate(sp, 3.7, 0)  == 3.231871225873364
        @test evaluate(sp, 3.7, 1)  == 0.5302076388029201
        @test evaluate(sp, 3.7, 2)  == -0.10845570338410908
        @test evaluate(sp, 4.55, 0) == 3.648200133827526
        @test evaluate(sp, 4.55, 1) == 0.4541838178776905
        @test evaluate(sp, 4.55, 2) == -0.07419003098143778
        @test evaluate(sp, 5.99, 0) == 4.237539062232049
        @test evaluate(sp, 5.99, 1) == 0.3712394302217441
        @test evaluate(sp, 5.99, 2) == -0.04496099505943021

        # Integrals
        @test evaluate_integral(sp, 1.0, 3.0) == 4.093271810090825
        @test evaluate_integral(sp, 2.0, 5.0) == 9.228837042885505
        @test evaluate_integral(sp, 1.5, 5.5) == 12.142267533433998
    end

    # ------------------------------------------------------------------
    # 2. Integer x input
    # ------------------------------------------------------------------
    @testset "Integer x spline" begin
        x = [1, 2, 3, 4, 5, 6, 7, 8]
        y = log.(x) .+ sqrt.(x)
        sp = Schumaker(x, y)

        @test evaluate(sp, 1.5) == 1.5808465653288954
        @test evaluate(sp, 3.5) == 3.122995892143333
        @test evaluate(sp, 6.5) == 4.42125103928065
    end

    # ------------------------------------------------------------------
    # 3. Extrapolation schemes
    # ------------------------------------------------------------------
    @testset "Extrapolation" begin
        x = collect(range(0.0, stop=10.0, length=40))
        y = x .^ 2

        sp_cc = Schumaker(x, y; extrapolation=(Curve, Curve))
        @test evaluate(sp_cc, -1.0) == 0.297089983173501
        @test evaluate(sp_cc,  0.0) == 0.0
        @test evaluate(sp_cc,  5.0) == 25.00003119419942
        @test evaluate(sp_cc, 10.0) == 100.0
        @test evaluate(sp_cc, 11.0) == 120.24019114435328

        sp_ll = Schumaker(x, y; extrapolation=(Linear, Linear))
        @test evaluate(sp_ll, -1.0) == -0.11539108708955387
        @test evaluate(sp_ll,  0.0) == 0.0
        @test evaluate(sp_ll, 11.0) == 119.8701123934658

        sp_const = Schumaker(x, y; extrapolation=(Constant, Constant))
        @test evaluate(sp_const, -1.0) == 0.0
        @test evaluate(sp_const, 11.0) == 100.0

        sp_cl = Schumaker(x, y; extrapolation=(Constant, Linear))
        @test evaluate(sp_cl, -1.0) == 0.0
        @test evaluate(sp_cl, 11.0) == 119.8701123934658
    end

    # ------------------------------------------------------------------
    # 4. Explicit gradients
    # ------------------------------------------------------------------
    @testset "Explicit gradients" begin
        x = collect(range(1.0, stop=6.0, length=100))
        y = log.(x) .+ sqrt.(x)
        grads = (1.0 ./ x) .+ (0.5 .* x .^ (-0.5))
        sp = Schumaker(x, y; gradients=grads)

        @test evaluate(sp, 1.5, 0) == 1.630209915891205
        @test evaluate(sp, 1.5, 1) == 1.0749370634210484
        @test evaluate(sp, 3.0, 0) == 2.8306630620851925
        @test evaluate(sp, 3.0, 1) == 0.6220052382669599
        @test evaluate(sp, 5.5, 0) == 4.049955973630384
        @test evaluate(sp, 5.5, 1) == 0.39501941067716073
    end

    # ------------------------------------------------------------------
    # 5. Roots and optima
    # ------------------------------------------------------------------
    @testset "Roots and optima" begin
        x = collect(range(0.0, stop=10.0, length=400))
        y = (x .- 3.0) .^ 2 .- 4.0
        sp = Schumaker(x, y)

        roots_info = find_roots(sp)
        @test length(roots_info.roots) == 2
        @test roots_info.roots[1]             == 1.0000000743222477
        @test roots_info.roots[2]             == 4.99999999600521
        @test roots_info.first_derivatives[1] == -4.000089134853394
        @test roots_info.first_derivatives[2] == 3.9998539244372373

        optima_info = find_optima(sp)
        @test length(optima_info.optima) == 1
        @test optima_info.optima[1] == 3.007518796992481
    end

    # ------------------------------------------------------------------
    # 6. Algebra operations
    # ------------------------------------------------------------------
    @testset "Algebra" begin
        x = collect(range(0.5, stop=10.0, length=40))
        y = x .^ 2
        sp = Schumaker(x, y)
        c = 7.6

        @test (sp + c)(1.0) == 8.600277282609294
        @test (sp * c)(1.0) == 7.602107347830637
        @test (sp / c)(1.0) == 0.1316154319222756
        @test (c - sp)(1.0) == 6.599722717390705

        @test (sp + c)(3.0) == 16.600316157738092
        @test (sp * c)(3.0) == 68.40240279880949
        @test (sp / c)(3.0) == 1.1842521260181702
        @test (c - sp)(3.0) == -1.4003161577380925

        @test (sp + c)(7.0) == 56.59988722853018
        @test (sp * c)(7.0) == 372.3991429368293
        @test (sp / c)(7.0) == 6.44735358270134
        @test (c - sp)(7.0) == -41.399887228530176
    end

    # ------------------------------------------------------------------
    # 7. Splice splines
    # ------------------------------------------------------------------
    @testset "Splice" begin
        x1 = collect(range(0.5, stop=10.0, length=40))
        y1 = x1 .^ 2
        x2 = [0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0]
        y2 = sqrt.(x2)
        sp1 = Schumaker(x1, y1)
        sp2 = Schumaker(x2, y2)
        spliced = splice_splines(sp1, sp2, 3.0)

        @test evaluate(spliced, 1.0) == 1.0002772826092945
        @test evaluate(spliced, 2.5) == 6.250361102364421
        @test evaluate(spliced, 3.0) == 1.7255674392384717
        @test evaluate(spliced, 5.0) == 2.2356790988309396
        @test evaluate(spliced, 8.0) == 2.8284271247461903
    end

    # ------------------------------------------------------------------
    # 8. 2D spline
    # ------------------------------------------------------------------
    @testset "2D spline" begin
        gridx = collect(1.0:10.0)
        gridy = collect(1.0:10.0)
        grid = gridx * gridy'
        sp2d = Schumaker2d(gridx, gridy, grid; extrapolation=(Linear, Linear))

        @test evaluate(sp2d, 3.0, 4.0) == 12.0
        @test evaluate(sp2d, 5.5, 5.5) == 30.25
        @test evaluate(sp2d, 2.3, 7.8) == 17.939999999999998
    end

    # ------------------------------------------------------------------
    # 9. Coefficient matrix snapshot (x^2 at integer points)
    # ------------------------------------------------------------------
    @testset "Coefficient matrix" begin
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [1.0, 4.0, 9.0, 16.0, 25.0]
        sp = Schumaker(x, y)

        @test sp.IntStarts_ == [1.0, 1.6666666666666667, 2.0, 2.6028416502090037,
                                3.0, 3.5726891826610885, 4.0, 4.333333333333333]
        @test sp.coefficient_matrix_ == [0.46291333609701096 2.3827822185373186 1.0;
                                         1.8516533443880445 3.0 2.7942607395124397;
                                         0.6349631257306276 4.234435562925363 4.0;
                                         1.4629436903335038 5.0 6.783451186751318;
                                         0.7316005100063907 6.162040603780009 9.0;
                                         1.31408720073865 7.0 12.768879137765417;
                                         1.3154289725932966 8.12304735160447 16.0;
                                         0.3288572431483237 9.0 18.85384122526741]
    end

    # ------------------------------------------------------------------
    # 10. Intersection points
    # ------------------------------------------------------------------
    @testset "Intersection points" begin
        x = collect(range(0.5, stop=15.0, length=60))
        y1 = x .^ 2
        y2 = 10.0 .* x
        s1 = Schumaker(x, y1)
        s2 = Schumaker(x, y2)
        ips = get_intersection_points(s1, s2)

        @test length(ips) == 1
        @test ips[1] == 10.00000752874307
    end

    # ------------------------------------------------------------------
    # 11. Derivative spline
    # ------------------------------------------------------------------
    @testset "Derivative spline" begin
        x = collect(range(1.0, stop=5.0, length=50))
        y = x .^ 3
        sp = Schumaker(x, y)
        dsp = find_derivative_spline(sp)

        @test evaluate(dsp, 1.5) == 6.764688163467375
        @test evaluate(dsp, 3.0) == 26.97405781968263
        @test evaluate(dsp, 4.5) == 60.763422546206094
    end

    # ------------------------------------------------------------------
    # 12. Numerical stability: no NaN when tsi ≈ x1 (alpha ≈ 0)
    # ------------------------------------------------------------------
    @testset "No NaN from near-zero alpha" begin
        # Construct gradients where g2 ≈ delta for some interval, pushing
        # tsi close to x1 and exercising the Inf*0 avoidance in
        # schumakerIndInterval!.
        x = [0.0, 1.0, 2.0, 3.0]
        y = [0.0, 1.0, 4.0, 9.0]
        # delta for interval [1,2] is 3.0. Setting g2=3.0 exactly makes
        # tsi = x1 for that interval in the Condition2 branch.
        grads = [0.5, 1.0, 3.0, 5.0]
        sp = Schumaker(x, y; gradients=grads)
        # The spline must be constructible and evaluate without NaN
        for pt in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
            val = evaluate(sp, pt)
            @test !isnan(val)
            @test !isinf(val)
        end
        # Also verify with gradients that are extremely close to delta
        # (floating point near-miss)
        delta = (y[3] - y[2]) / (x[3] - x[2])  # == 3.0
        grads2 = [0.5, 1.0, delta + 1e-15, 5.0]
        sp2 = Schumaker(x, y; gradients=grads2)
        for pt in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
            val = evaluate(sp2, pt)
            @test !isnan(val)
            @test !isinf(val)
        end
    end

    # ------------------------------------------------------------------
    # 13. reshape_values
    # ------------------------------------------------------------------
    @testset "reshape_values" begin
        xv = [1.0, 2.0, 3.0, 4.0, 5.0]
        yv = [1.0, 3.0, 2.5, 4.0, 3.8]
        reshaped = reshape_values(xv, yv; increasing=true, concave=true)

        @test reshaped == [1.0, 3.0, 3.0, 3.0, 3.0]
    end
end
