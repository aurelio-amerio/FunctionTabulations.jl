# FunctionTabulations

[![CI](https://github.com/aurelio-amerio/FunctionTabulations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aurelio-amerio/FunctionTabulations.jl/actions/workflows/CI.yml)
[![Coverage](https://codecov.io/gh/aurelio-amerio/FunctionTabulations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/aurelio-amerio/FunctionTabulations.jl)
[![][docs-stable-img]][docs-stable-url]

This package is a wrapper around [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl) to compute interpolation tables (FunctionTabulations) for functions of up to three variables with support for Unitful. 

This package stems from the common need in physics research projects to tabulate slow functions. With this package, it becomes possible to easily compute the values of a function and tabulate its values for up to three variables. 

# Installation
```julia
using Pkg; Pkg.add("FunctionTabulations")
```


# Usage: 1D FunctionTabulations
## Simple usage
To compute and load the tabulation of a function with one variable, it is necessary to use `create_tabulation_1D`:

```julia
using FunctionTabulations

func_1D(x) = sin(x)

sin_tabulation = create_tabulation_1D(func_1D, xmin = 0.0, xmax = 3.0, npoints = 100) # produces a file called `func_1D_data.jld2`

isapprox(sin(3.0), sin_tabulation(3.0), rtol=1e-3)
```
It is always necessary to specify the function to be tabulated, the interval and the number of points of the tabulation grid. 
The routine will compute `npoints` values of the function in the range `xmin:xmax`, save the table in a `.jld2` file, and return the interpolating function of the tabulated data. If a proper `.jld2` file already exists for that function, and it was computed in the same `x` range with the same number of points `npoints`, it will automatically be loaded on subsequent calls of `compute_tabulation_1D`.

## Logarithmic scale

Sometimes, it is more advantageous to tabulate a function using a logarithmic scale. For example, a function `y=f(x)` might display a more linear trend when the `x` or `f(x)` axis is expressed in logarithmic (`log10`) scale. In this case, it is possible to ask the routine to tabulate or interpolate that axis in logarithmic scale. For example, given the function `f(x)=10^x`, it is more advantageous to consider the `f(x)` axis in logarithmic scale for the interpolation:

```julia
func_1D(x) = 10^x

func_1D_tabulation = create_tabulation_1D(func_1D, xmin = 0.0, xmax = 3.0, npoints = 100, f_scale=:log10)

isapprox(func_1D(3.0), func_1D_tabulation(3.0), rtol=1e-3)
```

On the other hand, a function which needs to be interpolated in a wide range of `x` and does not vary harshly with `x` may benefit from the option `x_scale=:log10`.

## Unitful support 
FunctionTabulations.jl supports [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). It is possible to tabulate functions which have variables of the `Quantity` type, and/or return quantities. This means that units of measurement are fully supported.

```julia
using Unitful

func_1d(x) = x^2

itp_1d_1 = create_tabulation_1D(
    func_1d,
    custom_name = "1d_1",
    xmin = 0.0u"m",
    xmax = 3.0u"m",
    npoints = 100,
    x_scale = :linear,
    f_scale = :linear
)

isapprox(itp_1d_1(2.0u"m"), func_1d(2.0u"m"), rtol = 1e-3)
    
```
In the case of a tabulation with units of measurement, the unit of the result of the tabulated function is the one returned by `f(xmin)`. Nontheless, it supports automatic unit conversion for the input variables:

```julia
func_1d(x) = x^2

itp_1d_1 = create_tabulation_1D(
    func_1d,
    custom_name = "1d_1",
    xmin = 0.0u"m",
    xmax = 3.0u"m",
    npoints = 100,
    x_scale = :linear,
    f_scale = :linear
)

isapprox(itp_1d_1(2.0u"m"), itp_1d_1(200.0u"cm"), rtol = 1e-3)
    
```
## Save location
It is possible to specify the folder where to save/load the interpolation using the option `jld_base_path`. Furthermore it is possible to specify the name of the file in which the interpolation will be stored using the option `custom_name`. Please note that the suffix `_data.jld2` will always be appended to the filename when the tabulation is stored. 

```julia
func_1d(x) = sin(x)

itp_1d_2 = create_tabulation_1D(
    func_1d,
    jld_base_path = "interpolations",
    custom_name = "1d_2",
    xmin = 0.0,
    xmax = 3.0,
    npoints = 100,
) 
# will generate a file called "1d_2_data.jld2" in the folder "interpolations". 
# If the folder doesn't exist, the routine will create it.
```
## Interpolation options
It is possible to customize the behaviour of the interpolation of the tabulated grid using `interpolation_type` and `extrapolation_bc`. The default interpolation method is a linear interpolation (`:linear`), but it is also possible to use a cubic spline (`:cubic`). Furthermore, it is possible to specify the behaviour outside of the tabulation domain `(xmin, xmax)` using `extrapolation_bc`. The default behaviour is to throw an exception, but it is also possible for to extrapolate linearly:

```julia
using Interpolations: Line

func_1d(x) = sin(x)

itp_1d_1 = create_tabulation_1D(
    func_1d,
    xmin = 0.0,
    xmax = 3.0,
    npoints = 100,
    x_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Line,
)
```

For more information about the extrapolation conditions, see the [Interpolations.jl](https://juliamath.github.io/Interpolations.jl/latest/extrapolation/) documentation and pass the extrapolation function as the `extrapolation_bc` argument.

## Extra `args` and `kwargs`
In case it is necessary to pass further fixed positional arguments or keyword arguments to the tabulated function, `create_tabulation_1D` will pass the extra `args` and `kwargs` to the function to be tabulated:

```julia
func_1d(x, a; b) = x^2 + a*b

itp_1d_1 = create_tabulation_1D(
    func_1d,
    2, # the `a` positional argument 
    xmin = 0.0,
    xmax = 3.0,
    npoints = 100,
    b=1, # the `b` positional argument
)

isapprox(itp_1d_1(2.0), func_1d(2.0, 2; b=1), rtol = 1e-3)
```

# Usage: 2D FunctionTabulations & 3D FunctionTabulations
2D FunctionTabulations and 3D FunctionTabulations have a similar syntax, with extra parameters for the `y` and `z` variables:

```julia
func_2d(x, y) = sin(x) * sin(y)

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
isapprox(itp_2d_1(1.0, 1.3), func_2d(1.0, 1.3), rtol = 1e-3)
```

```julia
func_3d(x, y, z) = x * y + z

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

isapprox(itp_3d_1(1.0, 1.3, 2.5), func_3d(1.0, 1.3, 2.5), rtol = 1e-3)
```


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://aurelio-amerio.github.io/FunctionTabulations.jl