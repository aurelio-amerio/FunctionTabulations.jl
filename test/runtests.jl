using FunctionTabulations
using Unitful
using Test

@testset "create_tabulation_1D_no_units" begin
    func_1d(x) = sin(x)

    # create interpolations and test them

    itp_1d_1 = create_tabulation_1D(
        func_1d,
        xmin = 0.0,
        xmax = 3.0,
        npoints = 100,
        x_scale = :linear,
        f_scale = :linear
    )

    itp_1d_2 = create_tabulation_1D(
        func_1d,
        jld_base_path = "interpolations",
        custom_name = "1d_2",
        xmin = 1e-1,
        xmax = 3.0,
        npoints = 100,
        interpolation_type = :cubic,
        x_scale = :log10,
        f_scale = :log10
    )

    @test isapprox(itp_1d_1(2.0), func_1d(2.0), rtol = 1e-3)
    @test isapprox(itp_1d_2(2.0), func_1d(2.0), rtol = 1e-3)

    # load interpolations and test them

    itp_1d_1 = create_tabulation_1D(
        func_1d,
        xmin = 0.0,
        xmax = 3.0,
        npoints = 100,
        x_scale = :linear,
        f_scale = :linear
    )

    itp_1d_2 = create_tabulation_1D(
        func_1d,
        custom_name = "1d_2",
        jld_base_path = "interpolations",
        xmin = 1e-1,
        xmax = 3.0,
        npoints = 100,
        interpolation_type = :cubic,
        x_scale = :log10,
        f_scale = :log10
    )

    @test isapprox(itp_1d_1(2.0), func_1d(2.0), rtol = 1e-3)
    @test isapprox(itp_1d_2(2.0), func_1d(2.0), rtol = 1e-3)

    rm("func_1d_data.jld2")
    rm("interpolations/1d_2_data.jld2")
    rm("interpolations")
end

@testset "create_tabulation_1D_with_units" begin
    func_1d(x) = x^2

    # create interpolations and test them

    itp_1d_1 = create_tabulation_1D(
        func_1d,
        custom_name = "1d_1",
        xmin = 0.0u"m",
        xmax = 3.0u"m",
        npoints = 100,
        x_scale = :linear,
        f_scale = :linear
    )

    itp_1d_2 = create_tabulation_1D(
        func_1d,
        custom_name = "1d_2",
        xmin = 1e-1u"m",
        xmax = 3.0u"m",
        npoints = 100,
        x_scale = :log10,
        f_scale = :log10
    )

    @test isapprox(itp_1d_1(2.0u"m"), func_1d(2.0u"m"), rtol = 1e-3)
    @test isapprox(itp_1d_2(2.0u"m"), func_1d(2.0u"m"), rtol = 1e-3)

    # load interpolations and test them

    itp_1d_1 = create_tabulation_1D(
        func_1d,
        custom_name = "1d_1",
        xmin = 0.0u"m",
        xmax = 3.0u"m",
        npoints = 100,
        x_scale = :linear,
        f_scale = :linear
    )

    itp_1d_2 = create_tabulation_1D(
        func_1d,
        custom_name = "1d_2",
        xmin = 1e-1u"m",
        xmax = 3.0u"m",
        npoints = 100,
        x_scale = :log10,
        f_scale = :log10
    )

    @test isapprox(itp_1d_1(2.0u"m"), func_1d(2.0u"m"), rtol = 1e-3)
    @test isapprox(itp_1d_2(2.0u"m"), func_1d(2.0u"m"), rtol = 1e-3)

    rm("1d_1_data.jld2")
    rm("1d_2_data.jld2")
end

@testset "create_tabulation_2D_no_units" begin
    func_2d(x, y) = sin(x) * sin(y)

    # create interpolations and test them

    itp_2d_1 = create_tabulation_2D(
        func_2d,
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :linear,
        y_scale = :linear,
        f_scale = :linear
    )

    itp_2d_2 = create_tabulation_2D(
        func_2d,
        custom_name = "2d_2",
        jld_base_path = "interpolations",
        xmin = 1e-1,
        xmax = 1.0,
        ymin = 1e-1,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        interpolation_type = :cubic,
        x_scale = :log10,
        y_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_2d_1(1.0, 1.3), func_2d(1.0, 1.3), rtol = 1e-3)
    @test isapprox(itp_2d_2(1.0, 1.3), func_2d(1.0, 1.3), rtol = 1e-3)

    # load interpolations and test them

    itp_2d_1 = create_tabulation_2D(
        func_2d,
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :linear,
        y_scale = :linear,
        f_scale = :linear
    )

    itp_2d_2 = create_tabulation_2D(
        func_2d,
        custom_name = "2d_2",
        jld_base_path = "interpolations",
        xmin = 1e-1,
        xmax = 1.0,
        ymin = 1e-1,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        interpolation_type = :cubic,
        x_scale = :log10,
        y_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_2d_1(1.0, 1.3), func_2d(1.0, 1.3), rtol = 1e-3)
    @test isapprox(itp_2d_2(1.0, 1.3), func_2d(1.0, 1.3), rtol = 1e-3)

    rm("func_2d_data.jld2")
    rm("interpolations/2d_2_data.jld2")
    rm("interpolations")
end

@testset "create_tabulation_2D_with_units" begin
    func_2d(x, y) = x^2 / y

    # create interpolations and test them

    itp_2d_1 = create_tabulation_2D(
        func_2d,
        custom_name = "2d_1",
        xmin = 0.0u"m",
        xmax = 1.0u"m",
        ymin = 0.0u"s",
        ymax = 2.0u"s",
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :linear,
        y_scale = :linear,
        f_scale = :linear
    )

    itp_2d_2 = create_tabulation_2D(
        func_2d,
        custom_name = "2d_2",
        xmin = 1e-1u"m",
        xmax = 1.0u"m",
        ymin = 1e-1u"s",
        ymax = 2.0u"s",
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :log10,
        y_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_2d_1(1.0u"m", 1.3u"s"), func_2d(1.0u"m", 1.3u"s"), rtol = 1e-3)
    @test isapprox(itp_2d_2(1.0u"m", 1.3u"s"), func_2d(1.0u"m", 1.3u"s"), rtol = 1e-3)

    # load interpolations and test them

    itp_2d_1 = create_tabulation_2D(
        func_2d,
        custom_name = "2d_1",
        xmin = 0.0u"m",
        xmax = 1.0u"m",
        ymin = 0.0u"s",
        ymax = 2.0u"s",
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :linear,
        y_scale = :linear,
        f_scale = :linear
    )

    itp_2d_2 = create_tabulation_2D(
        func_2d,
        custom_name = "2d_2",
        xmin = 1e-1u"m",
        xmax = 1.0u"m",
        ymin = 1e-1u"s",
        ymax = 2.0u"s",
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :log10,
        y_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_2d_1(1.0u"m", 1.3u"s"), func_2d(1.0u"m", 1.3u"s"), rtol = 1e-3)
    @test isapprox(itp_2d_2(1.0u"m", 1.3u"s"), func_2d(1.0u"m", 1.3u"s"), rtol = 1e-3)

    rm("2d_1_data.jld2")
    rm("2d_2_data.jld2")
end

@testset "create_tabulation_3D_no_units" begin
    func_3d(x, y, z) = x * y + z

    # create interpolations and test them
    itp_3d_1 = create_tabulation_3D(
        func_3d,
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :linear,
        z_scale = :linear,
        f_scale = :linear
    )

    itp_3d_2 = create_tabulation_3D(
        func_3d,
        jld_base_path = "interpolations",
        custom_name = "3d_2",
        xmin = 1e-1,
        xmax = 1.0,
        ymin = 1e-1,
        ymax = 2.0,
        zmin = 1e-1,
        zmax = 3.0,
        interpolation_type = :cubic,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :log10,
        y_scale = :log10,
        z_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_3d_1(1.0, 1.3, 2.5), func_3d(1.0, 1.3, 2.5), rtol = 1e-3)
    @test isapprox(itp_3d_2(1.0, 1.3, 2.5), func_3d(1.0, 1.3, 2.5), rtol = 1e-2)

    # load interpolations and test them

    itp_3d_1 = create_tabulation_3D(
        func_3d,
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :linear,
        z_scale = :linear,
        f_scale = :linear
    )

    itp_3d_2 = create_tabulation_3D(
        func_3d,
        jld_base_path = "interpolations",
        custom_name = "3d_2",
        xmin = 1e-1,
        xmax = 1.0,
        ymin = 1e-1,
        ymax = 2.0,
        zmin = 1e-1,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        interpolation_type = :cubic,
        x_scale = :log10,
        y_scale = :log10,
        z_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_3d_1(1.0, 1.3, 2.5), func_3d(1.0, 1.3, 2.5), rtol = 1e-3)
    @test isapprox(itp_3d_2(1.0, 1.3, 2.5), func_3d(1.0, 1.3, 2.5), rtol = 1e-2)

    rm("func_3d_data.jld2")
    rm("interpolations/3d_2_data.jld2")
    rm("interpolations")
end


@testset "create_tabulation_3D_with_units" begin
    func_3d(x, y, z) = x * y + z

    # create interpolations and test them
    itp_3d_1 = create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0u"m",
        xmax = 1.0u"m",
        ymin = 0.0u"s^-1",
        ymax = 2.0u"s^-1",
        zmin = 0.0u"m/s",
        zmax = 3.0u"m/s",
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :linear,
        z_scale = :linear,
        f_scale = :linear
    )

    itp_3d_2 = create_tabulation_3D(
        func_3d,
        custom_name = "3d_2",
        xmin = 1e-1u"m",
        xmax = 1.0u"m",
        ymin = 1e-1u"s^-1",
        ymax = 2.0u"s^-1",
        zmin = 1e-1u"m/s",
        zmax = 3.0u"m/s",
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :log10,
        y_scale = :log10,
        z_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_3d_1(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-3)
    @test isapprox(itp_3d_2(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-2)

    # load interpolations and test them

    itp_3d_1 = create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0u"m",
        xmax = 1.0u"m",
        ymin = 0.0u"s^-1",
        ymax = 2.0u"s^-1",
        zmin = 0.0u"m/s",
        zmax = 3.0u"m/s",
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :linear,
        z_scale = :linear,
        f_scale = :linear
    )

    itp_3d_2 = create_tabulation_3D(
        func_3d,
        custom_name = "3d_2",
        xmin = 1e-1u"m",
        xmax = 1.0u"m",
        ymin = 1e-1u"s^-1",
        ymax = 2.0u"s^-1",
        zmin = 1e-1u"m/s",
        zmax = 3.0u"m/s",
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :log10,
        y_scale = :log10,
        z_scale = :log10,
        f_scale = :log10
    )
    @test isapprox(itp_3d_1(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-3)
    @test isapprox(itp_3d_2(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-2)

    rm("3d_1_data.jld2")
    rm("3d_2_data.jld2")
end

@testset "errors" begin
    func_1d(x) = sin(x)
    func_2d(x, y) = sin(x) * sin(y)
    func_3d(x, y, z) = x * y + z

    @test_throws ArgumentError FunctionTabulations.scaler(1.0, :nan)
    @test_throws ArgumentError FunctionTabulations.un_scaler(1.0, :nan)

    @test_throws ArgumentError create_tabulation_1D(
        func_1d,
        custom_name = "1d_1",
        xmin = 0.0,
        xmax = 3.0,
        npoints = 100,
        x_scale = :nan,
        f_scale = :linear
    )
    @test_throws ArgumentError create_tabulation_1D(
        func_1d,
        custom_name = "1d_1",
        xmin = 0.0,
        xmax = 3.0,
        npoints = 100,
        x_scale = :linear,
        f_scale = :nan
    )
    @test_throws ArgumentError create_tabulation_1D(
        func_1d,
        custom_name = "1d_1",
        xmin = 0.0,
        xmax = 3.0,
        npoints = 100,
        interpolation_type = :nan,
    )

    @test_throws ArgumentError create_tabulation_2D(
        func_2d,
        custom_name = "2d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :nan,
        y_scale = :linear,
        f_scale = :linear
    )
    @test_throws ArgumentError create_tabulation_2D(
        func_2d,
        custom_name = "2d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :linear,
        y_scale = :nan,
        f_scale = :linear
    )
    @test_throws ArgumentError create_tabulation_2D(
        func_2d,
        custom_name = "2d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        x_scale = :linear,
        y_scale = :linear,
        f_scale = :nan
    )
    @test_throws ArgumentError create_tabulation_2D(
        func_2d,
        custom_name = "2d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        npoints_x = 100,
        npoints_y = 200,
        interpolation_type = :nan,
    )

    @test_throws ArgumentError  create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :nan,
        y_scale = :linear,
        z_scale = :linear,
        f_scale = :linear
    )
    @test_throws ArgumentError  create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :nan,
        z_scale = :linear,
        f_scale = :linear
    )
    @test_throws ArgumentError  create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :linear,
        z_scale = :nan,
        f_scale = :linear
    )
    @test_throws ArgumentError  create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        x_scale = :linear,
        y_scale = :linear,
        z_scale = :linear,
        f_scale = :nan
    )
    @test_throws ArgumentError  create_tabulation_3D(
        func_3d,
        custom_name = "3d_1",
        xmin = 0.0,
        xmax = 1.0,
        ymin = 0.0,
        ymax = 2.0,
        zmin = 0.0,
        zmax = 3.0,
        npoints_x = 100,
        npoints_y = 200,
        npoints_z = 200,
        interpolation_type = :nan,
    )

    rm("1d_1_data.jld2")
    rm("2d_1_data.jld2")
    rm("3d_1_data.jld2")
end


