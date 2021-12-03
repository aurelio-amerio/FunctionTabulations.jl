using Pkg; Pkg.activate("./")
using Revise
using CachedInterpolations

using Interpolations
using JLD
using ProgressMeter
using Unitful
using Unitful: NoUnits
using Base.Threads


#%%
func_1d(x) = sin(x)

itp_1d = create_interpolation_1D(
    func_1d,
    xmin=0.0,
    xmax=3.0,
    npoints=100,
    scale_x=:linear,
    scale_f=:linear
)

isapprox(itp_1d(2.0),func_1d(2.0), rtol=1e-3)
#%%
func_2d(x,y) = sin(x)*cos(y)

itp_2d = create_interpolation_2D(
    func_2d,
    xmin=0.0,
    xmax=3.0,
    ymin=0.0,
    ymax=3.0,
    npoints_x=100,
    npoints_y=100,
    scale_x=:linear,
    scale_y=:linear,
    scale_f=:linear
)
isapprox(itp_2d(2.0, 1.3),func_2d(2.0, 1.3), rtol=1e-3)
#%%
itp_2d(2.0, 0.5)
func_2d(2.0, 1.3)