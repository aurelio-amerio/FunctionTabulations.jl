using Pkg;
Pkg.activate("./");
using Revise
using CachedInterpolations
using Unitful
using Unitful: FreeUnits
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
    xmin = 1e-1u"m",
    xmax = 3.0u"m",
    ymin = 1e-1u"s^-1",
    ymax = 3.0u"s^-1",
    zmin = 1e-1u"m/s",
    zmax = 3.0u"m/s",
    npoints_x = 500,
    npoints_y = 500,
    npoints_z = 500,
    scale_x = :log10,
    scale_y = :log10,
    scale_z = :log10,
    scale_f = :log10
)
isapprox(itp_3d_1(2.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(2.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-3)
isapprox(itp_3d_2(2.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(2.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-2)