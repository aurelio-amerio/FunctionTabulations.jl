using Pkg;
Pkg.activate("./");
using Revise
using CachedInterpolations
using Unitful
using Unitful: FreeUnits
using BenchmarkTools
#%%
func_3d(x, y, z) = x * y + z

# create interpolations and test them
itp_3d_1 = create_interpolation_3D(
    func_3d,
    custom_name = "3d_1",
    xmin = 0.0u"m",
    xmax = 3.0u"m",
    ymin = 0.0u"s^-1",
    ymax = 3.0u"s^-1",
    zmin = 0.0u"m/s",
    zmax = 3.0u"m/s",
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


@btime itp_3d_1(2.0u"m", 1.3u"s^-1", 2.5u"m/s")
@btime itp_3d_2(2.0, 1.3, 2.5)