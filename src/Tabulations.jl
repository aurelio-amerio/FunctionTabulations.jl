module Tabulations
using Interpolations
using JLD2
using ProgressMeter
using Unitful
using Unitful: FreeUnits
using Base.Threads

export create_tabulation_1D, create_tabulation_2D, create_tabulation_3D

# get_unit_annotation(x::Union{Quantity,FreeUnits}) = Quantity{T,dimension(x),W} where {T<:Real,W<:FreeUnits}


function scaler(x::Real, scale::Symbol)
    if scale == :linear
        return x
    elseif scale == :log10
        return log10(x)
    else
        throw(ArgumentError("Scale <$scale> is not supported"))
    end
end

function un_scaler(x::Real, scale::Symbol)
    if scale == :linear
        return x
    elseif scale == :log10
        return 10 .^ x
    else
        throw(ArgumentError("Scale <$scale> is not supported"))
    end
end

function create_tabulation_1D_no_units(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::Real,
    xmax::Real,
    npoints::Int,
    x_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
)

    func_name = nameof(func)

    arg(x) = func(x, args...; kwargs...)

    if isnothing(jld_base_path)
        base_path = pwd()
    else
        base_path = jld_base_path
    end

    if !(ispath(base_path))
        mkdir(base_path)
    end

    if isnothing(custom_name)
        filename = "$(func_name)_data.jld2"
    else
        filename = "$(custom_name)_data.jld2"
    end

    filepath = "$(base_path)/$(filename)"

    function load_file()
        if isfile(filepath)
            data = load(filepath)
            if data["xmin"] == xmin && data["xmax"] == xmax && data["npoints"] == npoints
                return true
            else
                return false
            end
        else
            return false
        end
    end

    if load_file()
        data = load(filepath)
        @info "$(filename) loaded!"
        x = data["x"]

        data_matrix = data["func"]

    else
        if x_scale == :linear
            x = range(xmin, xmax, length = npoints)
        elseif x_scale == :log10
            x = range(log10(xmin), log10(xmax), length = npoints)
        else
            throw(ArgumentError("X scale $x_scale not supported"))
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        p = Progress(Int(npoints))
        update!(p, 0)

        data_matrix = zeros(npoints)

        @threads for i = 1:npoints

            data_matrix[i] = arg(un_scaler(x[i], x_scale))

            next!(p)
        end

        data_dict = Dict{
            String,
            Union{
                String,
                Real,
                Array{Float64},
                StepRangeLen
            },
        }()
        data_dict["x"] = x
        data_dict["func"] = convert(Array{Float64}, data_matrix)

        data_dict["xmin"] = xmin
        data_dict["xmax"] = xmax
        data_dict["npoints"] = npoints

        save(filepath, data_dict)
        @info "$(filename) created and exported!"
    end

    if f_scale == :log10
        data_matrix[data_matrix.<=1e-300] .= 1e-300
    end

    knots = (x,)
    f_matrix = scaler.(data_matrix, f_scale)

    if interpolation_type == :linear
        itp = LinearInterpolation(knots, f_matrix, extrapolation_bc = extrapolation_bc())
    elseif interpolation_type == :cubic
        itp = CubicSplineInterpolation(knots, f_matrix, extrapolation_bc = extrapolation_bc())
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end


    function func_interp(x::Real)
        return un_scaler(itp(scaler(x, x_scale)), f_scale)::Float64
    end

    return func_interp::Function
end

function create_tabulation_2D_no_units(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::Real,
    xmax::Real,
    ymin::Real,
    ymax::Real,
    npoints_x::Int,
    npoints_y::Int,
    x_scale = :linear,
    y_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
)

    func_name = nameof(func)

    arg(x, y) = func(x, y, args...; kwargs...)

    if isnothing(jld_base_path)
        base_path = pwd()
    else
        base_path = jld_base_path
    end

    if !(ispath(base_path))
        mkdir(base_path)
    end

    if isnothing(custom_name)
        filename = "$(func_name)_data.jld2"
    else
        filename = "$(custom_name)_data.jld2"
    end

    filepath = "$(base_path)/$(filename)"

    function load_file()
        if isfile(filepath)
            data = load(filepath)
            if data["xmin"] == xmin && data["xmax"] == xmax && data["ymin"] == ymin && data["ymax"] == ymax && data["npoints_x"] == npoints_x && data["npoints_y"] == npoints_y
                return true
            else
                return false
            end
        else
            return false
        end
    end

    if load_file()
        data = load(filepath)
        @info "$(filename) loaded!"
        x = data["x"]
        y = data["y"]

        data_matrix = data["func"]

    else
        if x_scale == :linear
            x = range(xmin, xmax, length = npoints_x)
        elseif x_scale == :log10
            x = range(log10(xmin), log10(xmax), length = npoints_x)
        else
            throw(ArgumentError("X scale $x_scale not supported"))
        end

        if y_scale == :linear
            y = range(ymin, ymax, length = npoints_y)
        elseif y_scale == :log10
            y = range(log10(ymin), log10(ymax), length = npoints_y)
        else
            throw(ArgumentError("Y scale $y_scale not supported"))
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        p = Progress(Int(npoints_y))
        update!(p, 0)

        data_matrix = zeros(npoints_x, npoints_y)

        @threads for j = 1:npoints_y
            for i = 1:npoints_x
                data_matrix[i, j] = arg(un_scaler(x[i], x_scale), un_scaler(y[j], y_scale))
            end
            next!(p)
        end

        data_dict = Dict{
            String,
            Union{
                String,
                Real,
                Array{Float64},
                StepRangeLen,
            },
        }()
        data_dict["x"] = x
        data_dict["y"] = y
        data_dict["func"] = convert(Array{Float64}, data_matrix)

        data_dict["xmin"] = xmin
        data_dict["xmax"] = xmax
        data_dict["ymin"] = ymin
        data_dict["ymax"] = ymax
        data_dict["npoints_x"] = npoints_x
        data_dict["npoints_y"] = npoints_y

        save(filepath, data_dict)
        @info "$(filename) created and exported!"
    end

    if f_scale == :log10
        data_matrix[data_matrix.<=1e-300] .= 1e-300
    end

    knots = (x, y)
    f_matrix = scaler.(data_matrix, f_scale)

    if interpolation_type == :linear
        itp = LinearInterpolation(knots, f_matrix, extrapolation_bc = extrapolation_bc())
    elseif interpolation_type == :cubic
        itp = CubicSplineInterpolation(knots, f_matrix, extrapolation_bc = extrapolation_bc())
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end

    function func_interp(x::Real, y::Real)
        return un_scaler(itp(scaler(x, x_scale), scaler(y, y_scale)), f_scale)::Float64
    end

    return func_interp::Function
end


function create_tabulation_3D_no_units(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::Real,
    xmax::Real,
    ymin::Real,
    ymax::Real,
    zmin::Real,
    zmax::Real,
    npoints_x::Int,
    npoints_y::Int,
    npoints_z::Int,
    x_scale = :linear,
    y_scale = :linear,
    z_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
)

    func_name = nameof(func)

    arg(x, y, z) = func(x, y, z, args...; kwargs...)

    if isnothing(jld_base_path)
        base_path = pwd()
    else
        base_path = jld_base_path
    end

    if !(ispath(base_path))
        mkdir(base_path)
    end

    if isnothing(custom_name)
        filename = "$(func_name)_data.jld2"
    else
        filename = "$(custom_name)_data.jld2"
    end

    filepath = "$(base_path)/$(filename)"

    function load_file()
        if isfile(filepath)
            data = load(filepath)
            if data["xmin"] == xmin && data["xmax"] == xmax && data["ymin"] == ymin && data["ymax"] == ymax && data["zmin"] == zmin && data["zmax"] == zmax && data["npoints_x"] == npoints_x && data["npoints_y"] == npoints_y && data["npoints_z"] == npoints_z
                return true
            else
                return false
            end
        else
            return false
        end
    end

    if load_file()
        data = load(filepath)
        @info "$(filename) loaded!"
        x = data["x"]
        y = data["y"]
        z = data["z"]

        data_matrix = data["func"]

    else
        if x_scale == :linear
            x = range(xmin, xmax, length = npoints_x)
        elseif x_scale == :log10
            x = range(log10(xmin), log10(xmax), length = npoints_x)
        else
            throw(ArgumentError("X scale $x_scale not supported"))
        end

        if y_scale == :linear
            y = range(ymin, ymax, length = npoints_y)
        elseif y_scale == :log10
            y = range(log10(ymin), log10(ymax), length = npoints_y)
        else
            throw(ArgumentError("Y scale $y_scale not supported"))
        end

        if z_scale == :linear
            z = range(zmin, zmax, length = npoints_z)
        elseif z_scale == :log10
            z = range(log10(zmin), log10(zmax), length = npoints_z)
        else
            throw(ArgumentError("Z scale $z_scale not supported"))
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end

        p = Progress(Int(npoints_z * npoints_y))
        update!(p, 0)

        data_matrix = zeros(npoints_x, npoints_y, npoints_z)

        @threads for k = 1:npoints_z
            for j = 1:npoints_y
                for i = 1:npoints_x
                    data_matrix[i, j, k] = arg(un_scaler(x[i], x_scale), un_scaler(y[j], y_scale), un_scaler(z[k], z_scale))
                end
                next!(p)
            end
        end

        data_dict = Dict{
            String,
            Union{
                String,
                Real,
                Array{Float64},
                StepRangeLen
            },
        }()
        data_dict["x"] = x
        data_dict["y"] = y
        data_dict["z"] = z
        data_dict["func"] = convert(Array{Float64}, data_matrix)

        data_dict["xmin"] = xmin
        data_dict["xmax"] = xmax
        data_dict["ymin"] = ymin
        data_dict["ymax"] = ymax
        data_dict["zmin"] = zmin
        data_dict["zmax"] = zmax
        data_dict["npoints_x"] = npoints_x
        data_dict["npoints_y"] = npoints_y
        data_dict["npoints_z"] = npoints_z

        save(filepath, data_dict, compress = true)
        @info "$(filename) created and exported!"
    end

    if f_scale == :log10
        data_matrix[data_matrix.<=1e-300] .= 1e-300
    end

    knots = (x, y, z)
    f_matrix = scaler.(data_matrix, f_scale)

    if interpolation_type == :linear
        itp = LinearInterpolation(knots, f_matrix, extrapolation_bc = extrapolation_bc())
    elseif interpolation_type == :cubic
        itp = CubicSplineInterpolation(knots, f_matrix, extrapolation_bc = extrapolation_bc())
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end

    function func_interp(x::Real, y::Real, z::Real)
        return un_scaler(itp(scaler(x, x_scale), scaler(y, y_scale), scaler.(z, z_scale)), f_scale)::Float64
    end

    return func_interp::Function
end

function wrap_function_1D_add_units(func, x_units, f_units, args...; kwargs...)
    function wrapped_func(x)
        x_v = ustrip(x_units, x)
        return func(x_v, args...; kwargs...) * f_units
    end
    return wrapped_func
end

function wrap_function_2D_add_units(func, x_units, y_units, f_units, args...; kwargs...)
    function wrapped_func(x, y)
        x_v = ustrip(x_units, x)
        y_v = ustrip(y_units, y)
        return func(x_v, y_v, args...; kwargs...) * f_units
    end
    return wrapped_func
end

function wrap_function_3D_add_units(func, x_units, y_units, z_units, f_units, args...; kwargs...)
    function wrapped_func(x, y, z)
        x_v = ustrip(x_units, x)
        y_v = ustrip(y_units, y)
        z_v = ustrip(z_units, z)
        return func(x_v, y_v, z_v, args...; kwargs...) * f_units
    end
    return wrapped_func
end

function wrap_function_1D_remove_units(func, x_units, f_units, args...; kwargs...)
    function wrapped_function(x::Real)
        x_v = x * x_units
        return ustrip(f_units, func(x_v, args...; kwargs...))
    end
    return wrapped_function
end

function wrap_function_2D_remove_units(func, x_units, y_units, f_units, args...; kwargs...)
    function wrapped_function(x::Real, y::Real)
        x_v = x * x_units
        y_v = y * y_units
        return ustrip(f_units, func(x_v, y_v, args...; kwargs...))
    end
    return wrapped_function
end

function wrap_function_3D_remove_units(func, x_units, y_units, z_units, f_units, args...; kwargs...)
    function wrapped_function(x::Real, y::Real, z::Real)
        x_v = x * x_units
        y_v = y * y_units
        z_v = z * z_units
        return ustrip(f_units, func(x_v, y_v, z_v, args...; kwargs...))
    end
    return wrapped_function
end

##### interpolations with units ######

function create_tabulation_1D(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::T,
    xmax::T,
    npoints::Int,
    x_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Union{Real,Quantity}}

    x_units = unit(xmin)
    f_units = unit(func(xmin, args...; kwargs...))

    if f_units == NoUnits && x_units == NoUnits
        return create_tabulation_1D_no_units(
            func,
            args...;
            jld_base_path = jld_base_path,
            custom_name = custom_name,
            xmin = xmin,
            xmax = xmax,
            npoints = npoints,
            x_scale = x_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )
    else
        func_wrapped = wrap_function_1D_remove_units(func, x_units, f_units)

        xmin_v = ustrip(x_units, xmin)
        xmax_v = ustrip(x_units, xmax)

        if isnothing(custom_name)
            name = nameof(func)
        else
            name = custom_name
        end

        itp = create_tabulation_1D_no_units(
            func_wrapped,
            args...;
            jld_base_path = jld_base_path,
            custom_name = name,
            xmin = xmin_v,
            xmax = xmax_v,
            npoints = npoints,
            x_scale = x_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )

        scaled_itp = wrap_function_1D_add_units(itp, x_units, f_units)
        return scaled_itp
    end
end

function create_tabulation_2D(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::T,
    xmax::T,
    ymin::V,
    ymax::V,
    npoints_x::Int,
    npoints_y::Int,
    x_scale = :linear,
    y_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity}}

    x_units = unit(xmin)
    y_units = unit(ymin)
    f_units = unit(func(xmin, ymin, args...; kwargs...))

    if f_units == NoUnits && x_units == NoUnits && y_units == NoUnits
        return create_tabulation_2D_no_units(
            func,
            args...;
            jld_base_path = jld_base_path,
            custom_name = custom_name,
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            npoints_x = npoints_x,
            npoints_y = npoints_y,
            x_scale = x_scale,
            y_scale = y_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )
    else
        func_wrapped = wrap_function_2D_remove_units(func, x_units, y_units, f_units)

        xmin_v = ustrip(x_units, xmin)
        xmax_v = ustrip(x_units, xmax)
        ymin_v = ustrip(y_units, ymin)
        ymax_v = ustrip(y_units, ymax)

        if isnothing(custom_name)
            name = nameof(func)
        else
            name = custom_name
        end

        itp = create_tabulation_2D_no_units(
            func_wrapped,
            args...;
            jld_base_path = jld_base_path,
            custom_name = name,
            xmin = xmin_v,
            xmax = xmax_v,
            ymin = ymin_v,
            ymax = ymax_v,
            npoints_x = npoints_x,
            npoints_y = npoints_y,
            x_scale = x_scale,
            y_scale = y_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )

        scaled_itp = wrap_function_2D_add_units(itp, x_units, y_units, f_units)
        return scaled_itp
    end
end

function create_tabulation_3D(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::T,
    xmax::T,
    ymin::V,
    ymax::V,
    zmin::W,
    zmax::W,
    npoints_x::Int,
    npoints_y::Int,
    npoints_z::Int,
    x_scale = :linear,
    y_scale = :linear,
    z_scale = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity},W<:Union{Real,Quantity}}

    x_units = unit(xmin)
    y_units = unit(ymin)
    z_units = unit(zmin)
    f_units = unit(func(xmin, ymin, zmin, args...; kwargs...))

    if f_units == NoUnits && x_units == NoUnits && y_units == NoUnits && z_units == NoUnits
        return create_tabulation_3D_no_units(
            func,
            args...;
            jld_base_path = jld_base_path,
            custom_name = custom_name,
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            zmin = zmin,
            zmax = zmax,
            npoints_x = npoints_x,
            npoints_y = npoints_y,
            npoints_z = npoints_z,
            x_scale = x_scale,
            y_scale = y_scale,
            z_scale = z_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )
    else
        func_wrapped = wrap_function_3D_remove_units(func, x_units, y_units, z_units, f_units)

        xmin_v = ustrip(x_units, xmin)
        xmax_v = ustrip(x_units, xmax)
        ymin_v = ustrip(y_units, ymin)
        ymax_v = ustrip(y_units, ymax)
        zmin_v = ustrip(z_units, zmin)
        zmax_v = ustrip(z_units, zmax)

        if isnothing(custom_name)
            name = nameof(func)
        else
            name = custom_name
        end

        itp = create_tabulation_3D_no_units(
            func_wrapped,
            args...;
            jld_base_path = jld_base_path,
            custom_name = name,
            xmin = xmin_v,
            xmax = xmax_v,
            ymin = ymin_v,
            ymax = ymax_v,
            zmin = zmin_v,
            zmax = zmax_v,
            npoints_x = npoints_x,
            npoints_y = npoints_y,
            npoints_z = npoints_z,
            x_scale = x_scale,
            y_scale = y_scale,
            z_scale = z_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )

        scaled_itp = wrap_function_3D_add_units(itp, x_units, y_units, z_units, f_units)
        return scaled_itp
    end
end

end

# made with <3