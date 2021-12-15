module FunctionTabulations
using Interpolations
using JLD2
using ProgressMeter
using Unitful
using Unitful: FreeUnits
using Base.Threads

export create_tabulation_1D, create_tabulation_2D, create_tabulation_3D

# get_unit_annotation(x::Union{Quantity,FreeUnits}) = Quantity{T,dimension(x),W} where {T<:Real,W<:FreeUnits}

"""
    scaler(x::Real, scale::Symbol)

Utility function to scale a variable, accepts :linear and :log10 scale
"""
function scaler(x::Real, scale::Symbol)
    if scale == :linear
        return x
    elseif scale == :log10
        return log10(x)
    else
        throw(ArgumentError("Scale <$scale> is not supported"))
    end
end

"""
    un_scaler(x::Real, scale::Symbol)

Utility function to un_scale a variable, accepts :linear and :log10 scale
"""
function un_scaler(x::Real, scale::Symbol)
    if scale == :linear
        return x
    elseif scale == :log10
        return 10 .^ x
    else
        throw(ArgumentError("Scale <$scale> is not supported"))
    end
end

"""
Utility function to create an interpolation without units of measurement for 1D functions.
"""
function create_tabulation_1D_no_units(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::Real,
    xmax::Real,
    npoints::Int,
    x_scale::Symbol = :linear,
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

function create_tabulation_1D_no_units(
    func::Function,
    x::Vector{T},
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    x_scale::Symbol = :linear,
    f_scale = :linear,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Real}

    func_name = nameof(func)

    arg(x) = func(x, args...; kwargs...)

    x = scaler.(x, x_scale)
    npoints = length(x)
    xmin = minimum(x)
    xmax = maximum(x)

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
        throw(ArgumentError("$interpolation_type is not supported for this tabulation method"))
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end


    function func_interp(x::Real)
        return un_scaler(itp(scaler(x, x_scale)), f_scale)::Float64
    end

    return func_interp::Function
end

"""
Utility function to create an interpolation without units of measurement for 2D functions.
"""
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
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
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

function create_tabulation_2D_no_units(
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
    kwargs...
) where {T<:Real, V<:Real}

    func_name = nameof(func)

    arg(x, y) = func(x, y, args...; kwargs...)

    x = scaler.(x, x_scale)
    npoints_x = length(x)
    xmin = minimum(x)
    xmax = maximum(x)

    y = scaler.(y, y_scale)
    npoints_y = length(y)
    ymin = minimum(y)
    ymax = maximum(y)

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
        throw(ArgumentError("$interpolation_type is not supported for this tabulation method"))
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end

    function func_interp(x::Real, y::Real)
        return un_scaler(itp(scaler(x, x_scale), scaler(y, y_scale)), f_scale)::Float64
    end

    return func_interp::Function
end

"""
Utility function to create an interpolation without units of measurement for 3D functions.
"""
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
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    z_scale::Symbol = :linear,
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

function create_tabulation_3D_no_units(
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
    kwargs...
) where {T<:Real,V<:Real,W<:Real}

    func_name = nameof(func)

    arg(x, y, z) = func(x, y, z, args...; kwargs...)

    x = scaler.(x, x_scale)
    npoints_x = length(x)
    xmin = minimum(x)
    xmax = maximum(x)

    y = scaler.(y, y_scale)
    npoints_y = length(y)
    ymin = minimum(y)
    ymax = maximum(y)

    z = scaler.(z, z_scale)
    npoints_z = length(z)
    zmin = minimum(z)
    zmax = maximum(z)

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
        throw(ArgumentError("$interpolation_type is not supported for this tabulation method"))
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
- `jld_base_path`: path to the folder where the tabulation should be saved. Defaults to the current folder
- `custom_name`: custom name for the tabulation file, to which `_data.jld2` will be appended. Defaults to the name of the function to be tabulated
- `interpolation_type`: Type of the spline to be used for interpolation. Can either be :linear (1st order spline) or :cubic (3rd order spline)
- `extrapolation_bc`: behaviour of the tabulation outside of the boundaries defined by `[xmin, xmax]`. Possible behaviours are handled by Interpolations.jl and include `Throw`` (throws and eror if a value out of bounds is acessed) or `Line` (extrapolate linearly). Defaults to `Throw`
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
    xmin::T,
    xmax::T,
    npoints::Int,
    x_scale::Symbol = :linear,
    f_scale::Symbol = :linear,
    jld_base_path = nothing,
    custom_name = nothing,
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

function create_tabulation_1D(
    func::Function,
    x::Vector{T},
    args...;
    x_scale::Symbol = :linear,
    f_scale::Symbol = :linear,
    jld_base_path = nothing,
    custom_name = nothing,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Union{Real,Quantity}}
    xmin = minimum(x)

    x_units = unit(xmin)
    f_units = unit(func(xmin, args...; kwargs...))

    if f_units == NoUnits && x_units == NoUnits
        return create_tabulation_1D_no_units(
            func,
            x,
            args...;
            jld_base_path = jld_base_path,
            custom_name = custom_name,
            x_scale = x_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )
    else
        func_wrapped = wrap_function_1D_remove_units(func, x_units, f_units)

        x_v = ustrip.(x_units, x)

        if isnothing(custom_name)
            name = nameof(func)
        else
            name = custom_name
        end

        itp = create_tabulation_1D_no_units(
            func_wrapped,
            x_v,
            args...;
            jld_base_path = jld_base_path,
            custom_name = name,
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
    xmin::T,
    xmax::T,
    ymin::V,
    ymax::V,
    npoints_x::Int,
    npoints_y::Int,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    f_scale = :linear,
    jld_base_path = nothing,
    custom_name = nothing,
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

function create_tabulation_2D(
    func::Function,
    x::Vector{T},
    y::Vector{V},
    args...;
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    f_scale = :linear,
    jld_base_path = nothing,
    custom_name = nothing,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity}}
    xmin = minimum(x)
    ymin = minimum(y)
    x_units = unit(xmin)
    y_units = unit(ymin)
    f_units = unit(func(xmin, ymin, args...; kwargs...))

    if f_units == NoUnits && x_units == NoUnits && y_units == NoUnits
        return create_tabulation_2D_no_units(
            func,
            x,
            y,
            args...;
            jld_base_path = jld_base_path,
            custom_name = custom_name,
            x_scale = x_scale,
            y_scale = y_scale,
            f_scale = f_scale,
            interpolation_type = interpolation_type,
            extrapolation_bc = extrapolation_bc,
            kwargs...
        )
    else
        func_wrapped = wrap_function_2D_remove_units(func, x_units, y_units, f_units)

        x_v = ustrip.(x_units, x)
        y_v = ustrip.(y_units, y)

        if isnothing(custom_name)
            name = nameof(func)
        else
            name = custom_name
        end

        itp = create_tabulation_2D_no_units(
            func_wrapped,
            x_v,
            y_v,
            args...;
            jld_base_path = jld_base_path,
            custom_name = name,
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
    xmin::T,
    xmax::T,
    ymin::V,
    ymax::V,
    zmin::W,
    zmax::W,
    npoints_x::Int,
    npoints_y::Int,
    npoints_z::Int,
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    z_scale::Symbol = :linear,
    f_scale = :linear,
    jld_base_path = nothing,
    custom_name = nothing,
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

function create_tabulation_3D(
    func::Function,
    x::Vector{T},
    y::Vector{V},
    z::Vector{W},
    args...;
    x_scale::Symbol = :linear,
    y_scale::Symbol = :linear,
    z_scale::Symbol = :linear,
    f_scale = :linear,
    jld_base_path = nothing,
    custom_name = nothing,
    interpolation_type = :linear,
    extrapolation_bc = Throw,
    kwargs...
) where {T<:Union{Real,Quantity},V<:Union{Real,Quantity},W<:Union{Real,Quantity}}
    xmin = minimum(x)    
    ymin = minimum(y)
    zmin = minimum(z)

    x_units = unit(xmin)
    y_units = unit(ymin)
    z_units = unit(zmin)
    f_units = unit(func(xmin, ymin, zmin, args...; kwargs...))

    if f_units == NoUnits && x_units == NoUnits && y_units == NoUnits && z_units == NoUnits
        return create_tabulation_3D_no_units(
            func,
            x,
            y,
            z,
            args...;
            jld_base_path = jld_base_path,
            custom_name = custom_name,
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

        x_v = ustrip.(x_units, x)
        y_v = ustrip.(y_units, y)
        z_v = ustrip.(z_units, z)

        if isnothing(custom_name)
            name = nameof(func)
        else
            name = custom_name
        end

        itp = create_tabulation_3D_no_units(
            func_wrapped,
            x_v,
            y_v,
            z_v,
            args...;
            jld_base_path = jld_base_path,
            custom_name = name,
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