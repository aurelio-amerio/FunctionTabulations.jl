module FunctionTabulations
using Interpolations
using CodecZlib
using JLD2
using ProgressMeter
using Unitful
using Unitful: FreeUnits
using Base.Threads
using SHA: sha256


export create_tabulation_1D, create_tabulation_2D, create_tabulation_3D

# made with <3

include("helpers.jl")
include("helpers_1D.jl")
include("helpers_2D.jl")
include("helpers_3D.jl")


"""
    create_tabulation_1D(
        func::Function[,
        args...];
        xmin::T,
        xmax::T,
        npoints::Int[,
        x_scale::Symbol = :linear,
        f_scale = :linear,
        jld_base_path = nothing,
        custom_name = nothing,
        interpolation_type = :linear,
        extrapolation_bc = Throw,
        check_SHA = true,
        check_SHA_mode = :warn,
        metadata=nothing,
        metadata_validation_fn=nothing,
        kwargs...]
    ) where {T<:Union{Real,Quantity}}

Computes or loads the tabulation for a given function of the kind `f(x)`.

# Arguments
- `func::Function`: the function which should be tabulated. Must be of the kind `f(x)`
- `xmin::Union{Real,Quantity}`: the minimum `x` for the range which will be covered by the tabulation
- `xmax::Union{Real,Quantity}`: the maximum `x` for the range which will be covered by the tabulation
- `npoints::Int`: the number of points for which `f(x)` will be computed
- `x_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `x` points will be evenly spaced. If `:log10` is chosen, the points will be logarithmically spaced. Defaults to `:linear`
- `f_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `f(x)` points will not be alterated. If `:log10` is chosen, log10 will be applied to the `f(x)` points before interpolating. Defaults to `:linear` 
- `jld_base_path::String`: path to the folder where the tabulation should be saved. Defaults to the current folder
- `custom_name::String`: custom name for the tabulation file, to which `_data.jld2` will be appended. Defaults to the name of the function to be tabulated
- `interpolation_type::Symbol`: Type of the spline to be used for interpolation. Can either be :linear (1st order spline) or :cubic (3rd order spline)
- `extrapolation_bc`: behaviour of the tabulation outside of the boundaries defined by `[xmin, xmax]`. Possible behaviours are handled by Interpolations.jl and include `Throw`` (throws and eror if a value out of bounds is acessed) or `Line` (extrapolate linearly). Defaults to `Throw`
- `check_SHA::Bool`: Whether to check the SHA of the tabulated function. If the SHA of the tabulated function does not match the SHA stored inside the tabulation file, it might be necessary to recompute the tabulation, since it's likely that the definition of the tabulated function has changed.
- `check_SHA_mode::Symbol`: Describes the behaviour of this function if the SHA of the tabulated function does not match the SHA stored inside the tabulation file. Defaults to :warn (print a warning but load the tabulation all the same). Other options include :throw (throw an error to suggest manual deletion of the old tabulation) and :none (don't do anything).
- `metadata::Dict{String,Any}`: an optional `Dict{String, Any}` which contains additional metadata to be stored in the tabulation achive.
- `metadata_validation_fn::Function`: a function of the kind `metadata_validation_fn(metadata1::Dict{String, Any}, metadata2::Dict{String, Any})` to chech whether two metadata dictionaries contain the same data.
- `args..., kwargs...`: additional `args` and `kwargs` will be passed to the function which is to be tabulated

# Examples
```jldoctest
julia>  func_1d(x) = sin(x)

julia>  itp_1d_1 = create_tabulation_1D(
            func_1d,
            xmin = 0.0,
            xmax = 3.0,
            npoints = 100,
            x_scale = :linear,
            f_scale = :linear
        )

julia>  isapprox(itp_1d_1(2.0), func_1d(2.0), rtol = 1e-3)
true
```

Measurement units are supported too:

```jldoctest
julia>  using Unitful

julia>  func_1d(x) = x^2

julia>  itp_1d_1 = create_tabulation_1D(
            func_1d,
            xmin = 0.0u"m",
            xmax = 3.0u"m",
            npoints = 100,
            x_scale = :linear,
            f_scale = :linear
        )

julia>  isapprox(itp_1d_1(2.0u"m"), func_1d(2.0u"m"), rtol = 1e-3)
true
```
"""
function create_tabulation_1D(
    func::Function,
    args...;
    jld_base_path::Union{String,Nothing}=nothing,
    custom_name::Union{String,Nothing}=nothing,
    xmin::T,
    xmax::T,
    npoints::Int,
    x_scale::Symbol=:linear,
    f_scale::Symbol=:linear,
    interpolation_type::Symbol=:linear,
    extrapolation_bc=Throw,
    check_SHA::Bool=true,
    check_SHA_mode::Symbol=:warn,
    metadata::Union{Dict{String,Any},Nothing}=nothing,
    metadata_validation_fn::Union{Function,Nothing}=nothing,
    kwargs...
) where {T<:Union{Real,Quantity}}

    func_name = nameof(func)

    arg_tmp(x) = func(x, args...; kwargs...)
    x_units = unit(xmin)
    f_units = unit(arg_tmp(xmin))

    xmin_v = ustrip(x_units, xmin)::Real
    xmax_v = ustrip(x_units, xmax)::Real

    arg = wrap_function_1D_remove_units(arg_tmp, x_units, f_units)

    filename, filepath = _init(func_name, jld_base_path, custom_name)

    if _load_file_1D(filepath, xmin, xmax, npoints, metadata, metadata_validation_fn)
        data = load(filepath)

        @info "$(filename) loaded!"
        if check_SHA
            func_sha = data["sha"]
            test_sha(func, func_sha, mode=check_SHA_mode)
        end
        x = data["x"]

        data_matrix = data["func"]

    else
        if x_scale == :linear
            x = range(xmin_v, xmax_v, length=npoints)
        elseif x_scale == :log10
            x = range(log10(xmin_v), log10(xmax_v), length=npoints)
        else
            throw(ArgumentError("X scale $x_scale not supported"))
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        data_matrix = _compute_matrix_and_save_1D(func, arg, x, x_scale, xmin, xmax, npoints, metadata, filepath, filename)
    end

    func_interp_units = _get_interp_fn_1D(x, data_matrix, x_scale, f_scale, x_units, f_units, interpolation_type, extrapolation_bc)

    return func_interp_units::Function
end

"""
    create_tabulation_1D(
        func::Function,
        x::Vector{T},
        args...;
        jld_base_path = nothing,
        custom_name = nothing,
        x_scale::Symbol = :linear,
        f_scale = :linear,
        interpolation_type = :linear,
        extrapolation_bc = Throw,
        check_SHA = true,
        check_SHA_mode = :warn,
        metadata=nothing,
        metadata_validation_fn=nothing,
        kwargs...
        ) where {T<:Union{Real,Quantity}}

Provides an interface for `create_tabulation_1D` where the tabulation points are defined by the Vector `x`. Only `linear_interpolation` is supported.`
"""
function create_tabulation_1D(
    func::Function,
    x::Vector{T},
    args...;
    jld_base_path::Union{String,Nothing}=nothing,
    custom_name::Union{String,Nothing}=nothing,
    x_scale::Symbol=:linear,
    f_scale::Symbol=:linear,
    interpolation_type::Symbol=:linear,
    extrapolation_bc=Throw,
    check_SHA::Bool=true,
    check_SHA_mode::Symbol=:warn,
    metadata::Union{Dict{String,Any},Nothing}=nothing,
    metadata_validation_fn::Union{Function,Nothing}=nothing,
    kwargs...
) where {T<:Union{Real,Quantity}}

    func_name = nameof(func)

    arg_tmp(x) = func(x, args...; kwargs...)

    npoints = length(x)
    xmin = minimum(x)
    xmax = maximum(x)

    x_units = unit(xmin)
    f_units = unit(arg_tmp(xmin))

    x_v = ustrip.(x_units, x)
    x_v = scaler.(x_v, x_scale)

    arg = wrap_function_1D_remove_units(arg_tmp, x_units, f_units)

    filename, filepath = _init(func_name, jld_base_path, custom_name)

    if _load_file_1D(filepath, xmin, xmax, npoints, metadata, metadata_validation_fn)
        data = load(filepath)
        @info "$(filename) loaded!"
        if check_SHA
            func_sha = data["sha"]
            test_sha(func, func_sha, mode=check_SHA_mode)
        end
        x_v = data["x"]

        data_matrix = data["func"]

    else
        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        data_matrix = _compute_matrix_and_save_1D(func, arg, x_v, x_scale, xmin, xmax, npoints, metadata, filepath, filename)
    end


    func_interp_units = _get_interp_fn_1D(x_v, data_matrix, x_scale, f_scale, x_units, f_units, interpolation_type, extrapolation_bc)

    return func_interp_units::Function
end


"""
    create_tabulation_2D(
        func::Function[,
        args...];
        xmin::T,
        xmax::T,
        ymin::V,
        ymax::V,
        npoints_x::Int,
        npoints_y::Int[,
        x_scale::Symbol = :linear,
        y_scale::Symbol = :linear,
        f_scale = :linear,
        jld_base_path = nothing,
        custom_name = nothing,
        interpolation_type = :linear,
        extrapolation_bc = Throw,
        check_SHA = true,
        check_SHA_mode = :warn,
        metadata=nothing,
        metadata_validation_fn=nothing,
        kwargs...]
    ) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity}}

Computes or loads the tabulation for a given function of the kind `f(x, y)`.

# Arguments
- `func::Function`: the function which should be tabulated. Must be of the kind `f(x, y)`
- `xmin::Union{Real,Quantity}`: the minimum `x` for the range which will be covered by the tabulation
- `xmax::Union{Real,Quantity}`: the maximum `x` for the range which will be covered by the tabulation
- `ymin::Union{Real,Quantity}`: the minimum `y` for the range which will be covered by the tabulation
- `ymax::Union{Real,Quantity}`: the maximum `y` for the range which will be covered by the tabulation
- `npoints_x::Int`: the number of points along the `x` axis for which `f(x,y)` will be computed
- `npoints_y::Int`: the number of points along the `y` axis for which `f(x,y)` will be computed
- `x_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `x` points will be evenly spaced. If `:log10` is chosen, the points will be logarithmically spaced. Defaults to `:linear`
- `y_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `y` points will be evenly spaced. If `:log10` is chosen, the points will be logarithmically spaced. Defaults to `:linear`
- `f_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `f(x)` points will not be alterated. If `:log10` is chosen, log10 will be applied to the `f(x)` points before interpolating. Defaults to `:linear` 
- `jld_base_path`: path to the folder where the tabulation should be saved. Defaults to the current folder
- `custom_name`: custom name for the tabulation file, to which `_data.jld2` will be appended. Defaults to the name of the function to be tabulated
- `interpolation_type`: Type of the spline to be used for interpolation. Can either be :linear (1st order spline) or :cubic (3rd order spline)
- `extrapolation_bc`: behaviour of the tabulation outside of the boundaries defined by `[xmin, xmax]`. Possible behaviours are handled by Interpolations.jl and include `Throw`` (throws and eror if a value out of bounds is acessed) or `Line` (extrapolate linearly). Defaults to `Throw`
- `check_SHA::Bool`: Whether to check the SHA of the tabulated function. If the SHA of the tabulated function does not match the SHA stored inside the tabulation file, it might be necessary to recompute the tabulation, since it's likely that the definition of the tabulated function has changed.
- `check_SHA_mode::Symbol`: Describes the behaviour of this function if the SHA of the tabulated function does not match the SHA stored inside the tabulation file. Defaults to :warn (print a warning but load the tabulation all the same). Other options include :throw (throw an error to suggest manual deletion of the old tabulation) and :none (don't do anything).
- `metadata::Dict{String,Any}`: an optional `Dict{String, Any}` which contains additional metadata to be stored in the tabulation achive.
- `metadata_validation_fn::Function`: a function of the kind `metadata_validation_fn(metadata1::Dict{String, Any}, metadata2::Dict{String, Any})` to chech whether two metadata dictionaries contain the same data.
- `args..., kwargs...`: additional `args` and `kwargs` will be passed to the function which is to be tabulated

# Examples
```jldoctest
julia>  func_2d(x, y) = sin(x) * sin(y)

julia>  itp_2d_1 = create_tabulation_2D(
            func_2d,
            xmin = 0.0,
            xmax = 1.0,
            ymin = 0.0,
            ymax = 2.0,
            npoints_x = 100,
            npoints_y = 100,
        )

julia>  isapprox(itp_2d_1(1.0, 1.3), func_2d(1.0, 1.3), rtol = 1e-3)
true
```

Measurement units are supported too:

```jldoctest
julia>  using Unitful

julia>  func_2d(x, y) = x^2 + y

julia>  itp_2d_1 = create_tabulation_2D(
            func_2d,
            xmin = 0.0u"m",
            xmax = 1.0u"m",
            ymin = 0.0u"m^2",
            ymax = 2.0u"m^2",
            npoints_x = 100,
            npoints_y = 100,
        )

julia>  isapprox(itp_2d_1(1.0u"m", 1.3u"m^2"), func_2d(1.0u"m", 1.3u"m^2"), rtol = 1e-3)
true
```
"""
function create_tabulation_2D(
    func::Function,
    args...;
    jld_base_path::Union{String,Nothing}=nothing,
    custom_name::Union{String,Nothing}=nothing,
    xmin::T,
    xmax::T,
    ymin::V,
    ymax::V,
    npoints_x::Int,
    npoints_y::Int,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    f_scale::Symbol=:linear,
    interpolation_type::Symbol=:linear,
    extrapolation_bc=Throw,
    check_SHA::Bool=true,
    check_SHA_mode::Symbol=:warn,
    metadata::Union{Dict{String,Any},Nothing}=nothing,
    metadata_validation_fn::Union{Function,Nothing}=nothing,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity}}

    func_name = nameof(func)

    arg_tmp(x, y) = func(x, y, args...; kwargs...)
    x_units = unit(xmin)
    y_units = unit(ymin)
    f_units = unit(arg_tmp(xmin, ymin))

    xmin_v = ustrip(x_units, xmin)::Real
    xmax_v = ustrip(x_units, xmax)::Real
    ymin_v = ustrip(y_units, ymin)::Real
    ymax_v = ustrip(y_units, ymax)::Real

    arg = wrap_function_2D_remove_units(arg_tmp, x_units, y_units, f_units)

    filename, filepath = _init(func_name, jld_base_path, custom_name)

    if _load_file_2D(filepath, xmin, xmax, ymin, ymax, npoints_x, npoints_y, metadata, metadata_validation_fn)
        data = load(filepath)
        @info "$(filename) loaded!"
        if check_SHA
            func_sha = data["sha"]
            test_sha(func, func_sha, mode=check_SHA_mode)
        end
        x = data["x"]
        y = data["y"]

        data_matrix = data["func"]

    else
        if x_scale == :linear
            x = range(xmin_v, xmax_v, length=npoints_x)
        elseif x_scale == :log10
            x = range(log10(xmin_v), log10(xmax_v), length=npoints_x)
        else
            throw(ArgumentError("X scale $x_scale not supported"))
        end

        if y_scale == :linear
            y = range(ymin_v, ymax_v, length=npoints_y)
        elseif y_scale == :log10
            y = range(log10(ymin_v), log10(ymax_v), length=npoints_y)
        else
            throw(ArgumentError("Y scale $y_scale not supported"))
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        data_matrix = _compute_matrix_and_save_2D(func, arg, x, y, x_scale, y_scale, xmin, xmax, ymin, ymax, npoints_x, npoints_y, metadata, filepath, filename)
    end

    func_interp_units = _get_interp_fn_2D(x, y, data_matrix, x_scale, y_scale, f_scale, x_units, y_units, f_units, interpolation_type, extrapolation_bc)

    return func_interp_units::Function
end

"""
    function create_tabulation_2D(
        func::Function,
        x::Vector{T},
        y::Vector{V},
        args...;
        jld_base_path = nothing,
        custom_name = nothing,
        x_scale::Symbol = :linear,
        y_scale::Symbol = :linear,
        f_scale = :linear,
        interpolation_type = :linear,
        extrapolation_bc = Throw,
        check_SHA = true,
        check_SHA_mode = :warn,
        metadata=nothing,
        metadata_validation_fn=nothing,
        kwargs...
    ) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity}}

Provides an interface for `create_tabulation_2D` where the tabulation points are defined by the Vectors `x` and `y`. Only `linear_interpolation` is supported.`
"""
function create_tabulation_2D(
    func::Function,
    x::Vector{T},
    y::Vector{V},
    args...;
    jld_base_path::Union{String,Nothing}=nothing,
    custom_name::Union{String,Nothing}=nothing,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    f_scale::Symbol=:linear,
    interpolation_type::Symbol=:linear,
    extrapolation_bc=Throw,
    check_SHA::Bool=true,
    check_SHA_mode::Symbol=:warn,
    metadata::Union{Dict{String,Any},Nothing}=nothing,
    metadata_validation_fn::Union{Function,Nothing}=nothing,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity}}

    func_name = nameof(func)

    arg_tmp(x, y) = func(x, y, args...; kwargs...)

    npoints_x = length(x)
    xmin = minimum(x)
    xmax = maximum(x)

    npoints_y = length(y)
    ymin = minimum(y)
    ymax = maximum(y)

    x_units = unit(xmin)
    y_units = unit(ymin)
    f_units = unit(arg_tmp(xmin, ymin))

    x_v = ustrip.(x_units, x)
    x_v = scaler.(x_v, x_scale)
    y_v = ustrip.(y_units, y)
    y_v = scaler.(y_v, y_scale)

    arg = wrap_function_2D_remove_units(arg_tmp, x_units, y_units, f_units)

    filename, filepath = _init(func_name, jld_base_path, custom_name)

    if _load_file_2D(filepath, xmin, xmax, ymin, ymax, npoints_x, npoints_y, metadata, metadata_validation_fn)
        data = load(filepath)
        @info "$(filename) loaded!"
        if check_SHA
            func_sha = data["sha"]
            test_sha(func, func_sha, mode=check_SHA_mode)
        end
        x_v = data["x"]
        y_v = data["y"]

        data_matrix = data["func"]

    else
        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        data_matrix = _compute_matrix_and_save_2D(func, arg, x_v, y_v, x_scale, y_scale, xmin, xmax, ymin, ymax, npoints_x, npoints_y, metadata, filepath, filename)
    end

    func_interp_units = _get_interp_fn_2D(x_v, y_v, data_matrix, x_scale, y_scale, f_scale, x_units, y_units, f_units, interpolation_type, extrapolation_bc)

    return func_interp_units::Function
end


"""
    create_tabulation_3D(
        func::Function[,
        args...];
        xmin::T,
        xmax::T,
        ymin::V,
        ymax::V,
        zmin::W,
        zmax::W,
        npoints_x::Int,
        npoints_y::Int,
        npoints_z::Int[,
        x_scale::Symbol = :linear,
        y_scale::Symbol = :linear,
        z_scale::Symbol = :linear,
        f_scale = :linear,
        jld_base_path = nothing,
        custom_name = nothing,
        interpolation_type = :linear,
        extrapolation_bc = Throw,
        check_SHA = true,
        check_SHA_mode = :warn,
        metadata=nothing,
        metadata_validation_fn=nothing,
        kwargs...]
    ) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity},W<:Union{Real,Quantity}}

Computes or loads the tabulation for a given function of the kind `f(x, y, z)`.

# Arguments
- `func::Function`: the function which should be tabulated. Must be of the kind `f(x, y)`
- `xmin::Union{Real,Quantity}`: the minimum `x` for the range which will be covered by the tabulation
- `xmax::Union{Real,Quantity}`: the maximum `x` for the range which will be covered by the tabulation
- `ymin::Union{Real,Quantity}`: the minimum `y` for the range which will be covered by the tabulation
- `ymax::Union{Real,Quantity}`: the maximum `y` for the range which will be covered by the tabulation
- `zmin::Union{Real,Quantity}`: the minimum `z` for the range which will be covered by the tabulation
- `zmax::Union{Real,Quantity}`: the maximum `z` for the range which will be covered by the tabulation
- `npoints_x::Int`: the number of points along the `x` axis for which `f(x,y,z)` will be computed
- `npoints_y::Int`: the number of points along the `y` axis for which `f(x,y,z)` will be computed
- `npoints_z::Int`: the number of points along the `z` axis for which `f(x,y,z)` will be computed
- `x_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `x` points will be evenly spaced. If `:log10` is chosen, the points will be logarithmically spaced. Defaults to `:linear`
- `y_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `y` points will be evenly spaced. If `:log10` is chosen, the points will be logarithmically spaced. Defaults to `:linear`
- `z_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `z` points will be evenly spaced. If `:log10` is chosen, the points will be logarithmically spaced. Defaults to `:linear`
- `f_scale::Symbol`: the scale which will be used for tabulation. If `:linear` is choosen, the `f(x)` points will not be alterated. If `:log10` is chosen, log10 will be applied to the `f(x)` points before interpolating. Defaults to `:linear` 
- `jld_base_path`: path to the folder where the tabulation should be saved. Defaults to the current folder
- `custom_name`: custom name for the tabulation file, to which `_data.jld2` will be appended. Defaults to the name of the function to be tabulated
- `interpolation_type`: Type of the spline to be used for interpolation. Can either be :linear (1st order spline) or :cubic (3rd order spline)
- `extrapolation_bc`: behaviour of the tabulation outside of the boundaries defined by `[xmin, xmax]`. Possible behaviours are handled by Interpolations.jl and include `Throw`` (throws and eror if a value out of bounds is acessed) or `Line` (extrapolate linearly). Defaults to `Throw`
- `check_SHA::Bool`: Whether to check the SHA of the tabulated function. If the SHA of the tabulated function does not match the SHA stored inside the tabulation file, it might be necessary to recompute the tabulation, since it's likely that the definition of the tabulated function has changed.
- `check_SHA_mode::Symbol`: Describes the behaviour of this function if the SHA of the tabulated function does not match the SHA stored inside the tabulation file. Defaults to :warn (print a warning but load the tabulation all the same). Other options include :throw (throw an error to suggest manual deletion of the old tabulation) and :none (don't do anything).
- `metadata::Dict{String,Any}`: an optional `Dict{String, Any}` which contains additional metadata to be stored in the tabulation achive.
- `metadata_validation_fn::Function`: a function of the kind `metadata_validation_fn(metadata1::Dict{String, Any}, metadata2::Dict{String, Any})` to chech whether two metadata dictionaries contain the same data.
- `args..., kwargs...`: additional `args` and `kwargs` will be passed to the function which is to be tabulated

# Examples
```jldoctest
julia>  func_3d(x, y, z) = x * y + z

julia>  itp_3d_1 = create_tabulation_3D(
            func_3d,
            xmin = 0.0,
            xmax = 1.0,
            ymin = 0.0,
            ymax = 2.0,
            zmin = 0.0,
            zmax = 3.0,
            npoints_x = 100,
            npoints_y = 100,
            npoints_z = 100,
        )

julia>  isapprox(itp_3d_1(1.0, 1.3, 2.5), func_3d(1.0, 1.3, 2.5), rtol = 1e-3)
true
```

Measurement units are supported too:

```jldoctest
julia>  using Unitful

julia>  func_3d(x, y, z) = x * y + z

julia>  itp_3d_1 = create_tabulation_3D(
            func_3d,
            xmin = 0.0u"m",
            xmax = 1.0u"m",
            ymin = 0.0"s^-1",
            ymax = 2.0"s^-1",
            zmin = 0.0u"m/s",
            zmax = 3.0u"m/s",
            npoints_x = 100,
            npoints_y = 100,
            npoints_z = 100,
        )

julia>  isapprox(itp_3d_1(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), func_3d(1.0u"m", 1.3u"s^-1", 2.5u"m/s"), rtol = 1e-3)
true
```
"""
function create_tabulation_3D(
    func::Function,
    args...;
    jld_base_path::Union{String,Nothing}=nothing,
    custom_name::Union{String,Nothing}=nothing,
    xmin::T,
    xmax::T,
    ymin::V,
    ymax::V,
    zmin::W,
    zmax::W,
    npoints_x::Int,
    npoints_y::Int,
    npoints_z::Int,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    z_scale::Symbol=:linear,
    f_scale::Symbol=:linear,
    interpolation_type::Symbol=:linear,
    extrapolation_bc=Throw,
    check_SHA::Bool=true,
    check_SHA_mode::Symbol=:warn,
    compress::Bool=false,
    metadata::Union{Dict{String,Any},Nothing}=nothing,
    metadata_validation_fn::Union{Function,Nothing}=nothing,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity},W<:Union{Real,Quantity}}

    func_name = nameof(func)

    arg_tmp(x, y, z) = func(x, y, z, args...; kwargs...)

    x_units = unit(xmin)
    y_units = unit(ymin)
    z_units = unit(zmin)
    f_units = unit(arg_tmp(xmin, ymin, zmin))

    xmin_v = ustrip(x_units, xmin)::Real
    xmax_v = ustrip(x_units, xmax)::Real
    ymin_v = ustrip(y_units, ymin)::Real
    ymax_v = ustrip(y_units, ymax)::Real
    zmin_v = ustrip(z_units, zmin)::Real
    zmax_v = ustrip(z_units, zmax)::Real

    arg = wrap_function_3D_remove_units(arg_tmp, x_units, y_units, z_units, f_units)

    filename, filepath = _init(func_name, jld_base_path, custom_name)


    if _load_file_3D(filepath, xmin, xmax, ymin, ymax, zmin, zmax, npoints_x, npoints_y, npoints_z, metadata, metadata_validation_fn)
        data = load(filepath)
        @info "$(filename) loaded!"
        if check_SHA
            func_sha = data["sha"]
            test_sha(func, func_sha, mode=check_SHA_mode)
        end
        x = data["x"]
        y = data["y"]
        z = data["z"]

        data_matrix = data["func"]

    else
        if x_scale == :linear
            x = range(xmin_v, xmax_v, length=npoints_x)
        elseif x_scale == :log10
            x = range(log10(xmin_v), log10(xmax_v), length=npoints_x)
        else
            throw(ArgumentError("X scale $x_scale not supported"))
        end

        if y_scale == :linear
            y = range(ymin_v, ymax_v, length=npoints_y)
        elseif y_scale == :log10
            y = range(log10(ymin_v), log10(ymax_v), length=npoints_y)
        else
            throw(ArgumentError("Y scale $y_scale not supported"))
        end

        if z_scale == :linear
            z = range(zmin_v, zmax_v, length=npoints_z)
        elseif z_scale == :log10
            z = range(log10(zmin_v), log10(zmax_v), length=npoints_z)
        else
            throw(ArgumentError("Z scale $z_scale not supported"))
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        data_matrix = _compute_matrix_and_save_3D(func, arg, x, y, z, x_scale, y_scale, z_scale, xmin, xmax, ymin, ymax, zmin, zmax, npoints_x, npoints_y, npoints_z, metadata, filepath, filename, compress)
    end

    func_interp_units = _get_interp_fn_3D(x, y, z, data_matrix, x_scale, y_scale, z_scale, f_scale, x_units, y_units, z_units, f_units, interpolation_type, extrapolation_bc)

    return func_interp_units::Function
end

"""
    function create_tabulation_3D(
        func::Function,
        x::Vector{T},
        y::Vector{V},
        z::Vector{W},
        args...;
        jld_base_path = nothing,
        custom_name = nothing,
        x_scale::Symbol = :linear,
        y_scale::Symbol = :linear,
        z_scale::Symbol = :linear,
        f_scale = :linear,
        interpolation_type = :linear,
        extrapolation_bc = Throw,
        check_SHA = true,
        check_SHA_mode = :warn,
        metadata=nothing,
        metadata_validation_fn=nothing,
        kwargs...
    ) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity},W<:Union{Real,Quantity}}

Provides an interface for `create_tabulation_3D` where the tabulation points are defined by the Vectors `x`, `y` and `z`. Only `linear_interpolation` is supported.`
"""
function create_tabulation_3D(
    func::Function,
    x::Vector{T},
    y::Vector{V},
    z::Vector{W},
    args...;
    jld_base_path=nothing,
    custom_name=nothing,
    x_scale::Symbol=:linear,
    y_scale::Symbol=:linear,
    z_scale::Symbol=:linear,
    f_scale=:linear,
    interpolation_type=:linear,
    extrapolation_bc=Throw,
    check_SHA::Bool=true,
    check_SHA_mode::Symbol=:warn,
    compress::Bool=false,
    metadata::Union{Dict{String,Any},Nothing}=nothing,
    metadata_validation_fn::Union{Function,Nothing}=nothing,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity},W<:Union{Real,Quantity}}

    func_name = nameof(func)

    arg_tmp(x, y, z) = func(x, y, z, args...; kwargs...)

    npoints_x = length(x)
    xmin = minimum(x)
    xmax = maximum(x)

    npoints_y = length(y)
    ymin = minimum(y)
    ymax = maximum(y)

    npoints_z = length(z)
    zmin = minimum(z)
    zmax = maximum(z)

    x_units = unit(xmin)
    y_units = unit(ymin)
    z_units = unit(zmin)
    f_units = unit(arg_tmp(xmin, ymin, zmin))

    x_v = ustrip.(x_units, x)
    x_v = scaler.(x_v, x_scale)
    y_v = ustrip.(y_units, y)
    y_v = scaler.(y_v, y_scale)
    z_v = ustrip.(z_units, z)
    z_v = scaler.(z_v, z_scale)

    arg = wrap_function_3D_remove_units(arg_tmp, x_units, y_units, z_units, f_units)


    filename, filepath = _init(func_name, jld_base_path, custom_name)

    if _load_file_3D(filepath, xmin, xmax, ymin, ymax, zmin, zmax, npoints_x, npoints_y, npoints_z, metadata, metadata_validation_fn)
        data = load(filepath)
        @info "$(filename) loaded!"
        if check_SHA
            func_sha = data["sha"]
            test_sha(func, func_sha, mode=check_SHA_mode)
        end
        x_v = data["x"]
        y_v = data["y"]
        z_v = data["z"]

        data_matrix = data["func"]

    else

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        data_matrix = _compute_matrix_and_save_3D(func, arg, x_v, y_v, z_v, x_scale, y_scale, z_scale, xmin, xmax, ymin, ymax, zmin, zmax, npoints_x, npoints_y, npoints_z, metadata, filepath, filename, compress)
    end

    func_interp_units = _get_interp_fn_3D(x_v, y_v, z_v, data_matrix, x_scale, y_scale, z_scale, f_scale, x_units, y_units, z_units, f_units, interpolation_type, extrapolation_bc)

    return func_interp_units::Function
end

end