using Pkg;
Pkg.activate("./");
using Revise
using CachedInterpolations
#%%


func_1d(x) = sin(x)

itp_1d_1 = create_interpolation_1D(
    func_1d,
    custom_name = "1d_1",
    xmin = 0.0,
    xmax = 3.0,
    npoints = 100,
    scale_x = :linear,
    scale_f = :linear
)

itp_1d_2 = create_interpolation_1D(
    func_1d,
    custom_name = "1d_2",
    xmin = 1e-1,
    xmax = 3.0,
    npoints = 100,
    scale_x = :log10,
    scale_f = :log10
)

isapprox(itp_1d_1(2.0), func_1d(2.0), rtol = 1e-3)
isapprox(itp_1d_2(2.0), func_1d(2.0), rtol = 1e-3)
#%%



func_2d(x, y) = sin(x) * sin(y)

itp_2d_1 = create_interpolation_2D(
    func_2d,
    custom_name = "2d_1",
    xmin = 0.0,
    xmax = 3.0,
    ymin = 0.0,
    ymax = 3.0,
    npoints_x = 100,
    npoints_y = 200,
    scale_x = :linear,
    scale_y = :linear,
    scale_f = :linear
)

itp_2d_2 = create_interpolation_2D(
    func_2d,
    custom_name = "2d_2",
    xmin = 1e-1,
    xmax = 3.0,
    ymin = 1e-1,
    ymax = 3.0,
    npoints_x = 100,
    npoints_y = 200,
    scale_x = :log10,
    scale_y = :log10,
    scale_f = :log10
)
isapprox(itp_2d_1(2.0, 1.3), func_2d(2.0, 1.3), rtol = 1e-3)
isapprox(itp_2d_2(2.0, 1.3), func_2d(2.0, 1.3), rtol = 1e-3)


func_3d(x, y, z) = sin(x) * sin(y) + cos(z)

itp_3d_1 = create_interpolation_3D(
    func_3d,
    custom_name = "3d_1",
    xmin = 0.0,
    xmax = 3.0,
    ymin = 0.0,
    ymax = 3.0,
    zmin = 0.0,
    zmax = 3.0,
    npoints_x = 100,
    npoints_y = 200,
    npoints_z = 200,
    scale_x = :linear,
    scale_y = :linear,
    scale_z = :linear,
    scale_f = :linear
)

itp_3d_2 = create_interpolation_3D(
    func_3d,
    custom_name = "3d_2",
    xmin = 1e-1,
    xmax = 3.0,
    ymin = 1e-1,
    ymax = 3.0,
    zmin = 1e-1,
    zmax = 3.0,
    npoints_x = 500,
    npoints_y = 500,
    npoints_z = 500,
    scale_x = :log10,
    scale_y = :log10,
    scale_z = :log10,
    scale_f = :log10
)
isapprox(itp_3d_1(2.0, 1.3, 2.5), func_3d(2.0, 1.3, 2.5), rtol = 1e-3)
isapprox(itp_3d_2(2.0, 1.3, 2.5), func_3d(2.0, 1.3, 2.5), rtol = 1e-2)

itp_3d_2(2.0, 1.3, 2.5)
func_3d(2.0, 1.3, 2.5)