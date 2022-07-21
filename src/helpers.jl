
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
compute_SHA(func)

Utility function to compute the SHA of a function.
"""
function compute_SHA(func::Function)
    c1 = code_lowered(func)
    h1 = bytes2hex(sha256(string(c1)))
    return h1
end

"""
test_sha(func, sha[; mode = :warn])

Compares the SHA of the `func` funtion with the provided SHA.
"""
function test_sha(func::Function, sha::String; mode::Symbol=:warn)
    func_sha = compute_SHA(func)
    func_name = nameof(func)
    if sha != func_sha
        if mode == :warn
            @warn "The SHA for `$func_name` did not match the one of the stored tabulated function. Please check if the function definition has changed."
        elseif mode == :throw
            @error "The SHA for `$func_name` did not match the one of the stored tabulated function. Please check if the function definition has changed."
        elseif mode == :none
        else
            @error "Mode not supported. Got $(mode)."
        end
    end
end

"""
wrap_function_1D_add_units(func, x_units, f_units[, args...; kwargs...])

Add units back to a 1D function.
"""
function wrap_function_1D_add_units(func::Function, x_units, f_units, args...; kwargs...)
    function wrapped_func(x)
        x_v = ustrip(x_units, x)
        return func(x_v, args...; kwargs...) * f_units
    end
    return wrapped_func
end

"""
wrap_function_2D_add_units(func, x_units, y_units, f_units[, args...; kwargs...])

Add units back to a 2D function.
"""
function wrap_function_2D_add_units(func::Function, x_units, y_units, f_units, args...; kwargs...)
    function wrapped_func(x, y)
        x_v = ustrip(x_units, x)
        y_v = ustrip(y_units, y)
        return func(x_v, y_v, args...; kwargs...) * f_units
    end
    return wrapped_func
end

"""
wrap_function_3D_add_units(func, x_units, y_units, z_units, f_units[, args...; kwargs...])

Add units back to a 3D function.
"""
function wrap_function_3D_add_units(func::Function, x_units, y_units, z_units, f_units, args...; kwargs...)
    function wrapped_func(x, y, z)
        x_v = ustrip(x_units, x)
        y_v = ustrip(y_units, y)
        z_v = ustrip(z_units, z)
        return func(x_v, y_v, z_v, args...; kwargs...) * f_units
    end
    return wrapped_func
end

"""
wrap_function_1D_remove_units(func, x_units, f_units[, args...; kwargs...])

Remove units from a 1D function.
"""
function wrap_function_1D_remove_units(func::Function, x_units, f_units, args...; kwargs...)
    function wrapped_function(x::Real)
        x_v = x * x_units
        return ustrip(f_units, func(x_v, args...; kwargs...))
    end
    return wrapped_function
end

"""
wrap_function_2D_remove_units(func, x_units, y_units, f_units[, args...; kwargs...])

Remove units from a 2D function.
"""
function wrap_function_2D_remove_units(func::Function, x_units, y_units, f_units, args...; kwargs...)
    function wrapped_function(x::Real, y::Real)
        x_v = x * x_units
        y_v = y * y_units
        return ustrip(f_units, func(x_v, y_v, args...; kwargs...))
    end
    return wrapped_function
end

"""
wrap_function_3D_remove_units(func, x_units, y_units, z_units, f_units[, args...; kwargs...])

Remove units from a 3D function.
"""
function wrap_function_3D_remove_units(func::Function, x_units, y_units, z_units, f_units, args...; kwargs...)
    function wrapped_function(x::Real, y::Real, z::Real)
        x_v = x * x_units
        y_v = y * y_units
        z_v = z * z_units
        return ustrip(f_units, func(x_v, y_v, z_v, args...; kwargs...))
    end
    return wrapped_function
end

function _init(func_name, jld_base_path, custom_name)
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

    return filename, filepath
end


function get_metadata(filepath)
    data = load(filepath)
    return data["metadata"]
end

function check_metadata(metadata, loaded_metadata, metadata_validation_fn::Union{Function, Nothing}=nothing)   
    if metadata_validation_fn === nothing
        return metadata == loaded_metadata
    else
        return metadata_validation_fn(metadata, loaded_metadata)
    end
end
    
