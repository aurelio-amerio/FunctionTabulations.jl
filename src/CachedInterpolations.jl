module CachedInterpolations
using Interpolations
using JLD
using ProgressMeter

using Base.Threads

export create_interpolation_1D, create_interpolation_2D, create_interpolation_3D

function create_interpolation_1D(
    func::Function,
    args...;
    jld_base_path = nothing,
    custom_name = nothing,
    xmin::Real,
    xmax::Real,
    npoints::Int,
    scale_x=:linear,
    scale_f=:linear,
    extrapolation_bc = Throw,
    kwargs...,
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
        filename = "$(func_name)_data.jld"
    else
        filename = "$(custom_name)_data.jld"
    end
    
    filepath = "$(base_path)/$(filename)"

    function load_file()
        if isfile(filepath)
            data = load(filepath)
            if data["xmin"] == xmin && data["xmax"] == xmax 
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
        if scale_x == :linear
            x = range(xmin, xmax, length = npoints)
        elseif scale_x == :log
            x = 10 .^ range(log10(xmin), log10(xmax), length = npoints)
        else
            error("X scale $scale_x not supported")
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

            data_matrix[i] = arg(x[i])

            next!(p)
        end

        data_dict = Dict{
            String,
            Union{
                String,
                Real,
                Array{Float64},
            },
        }()
        data_dict["x"] = collect(x)
        data_dict["func"] = convert(Array{Float64}, data_matrix)

        data_dict["xmin"] = xmin
        data_dict["xmax"] = xmax

        save(filepath, data_dict)
        @info "$(filename) created and exported!"
    end

    data_matrix[data_matrix.<1e-299] .= 1e-300

    function scaler(x::Real, scale::Symbol)
        if scale == :linear
            return x
        elseif scale == :log
            return log10(x)
        else
            error("Scale <$scale> is not supported")
        end
    end

    function un_scaler(x::Real, scale::Symbol)
        if scale == :linear
            return x
        elseif scale == :log
            return 10 .^ x
        else
            error("Scale <$scale> is not supported")
        end
    end


    knots = (scaler.(x, scale_x), )
    f_matrix = scaler.(data_matrix, scale_f)

    itp = LinearInterpolation(knots, f_matrix , extrapolation_bc = extrapolation_bc())

    function func_interp(x::Real)
        return un_scaler(itp(scaler(x, scale_x)), scale_f)::Float64
    end

    return func_interp::Function
end

function create_interpolation_2D(
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
    scale_x=:linear,
    scale_y=:linear,
    scale_f=:linear,
    extrapolation_bc = Throw,
    kwargs...,
)

    func_name = nameof(func)

    arg(x,y) = func(x, y, args...; kwargs...)

    if isnothing(jld_base_path)  
        base_path = pwd()
    else
        base_path = jld_base_path
    end

    if !(ispath(base_path))
        mkdir(base_path)
    end

    if isnothing(custom_name)
        filename = "$(func_name)_data.jld"
    else
        filename = "$(custom_name)_data.jld"
    end
    
    filepath = "$(base_path)/$(filename)"

    function load_file()
        if isfile(filepath)
            data = load(filepath)
            if data["xmin"] == xmin && data["xmax"] == xmax && data["ymin"] == ymin && data["ymax"] == ymax 
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
        if scale_x == :linear
            x = range(xmin, xmax, length = npoints_x)
        elseif scale_x == :log
            x = 10 .^ range(log10(xmin), log10(xmax), length = npoints_x)
        else
            error("X scale $scale_x not supported")
        end

        if scale_y == :linear
            y = range(ymin, ymax, length = npoints_y)
        elseif scale_y == :log
            y = 10 .^ range(log10(ymin), log10(ymax), length = npoints_y)
        else
            error("Y scale $scale_y not supported")
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
                data_matrix[i] = arg(x[i], y[j])
            end
            next!(p)
        end

        data_dict = Dict{
            String,
            Union{
                String,
                Real,
                Array{Float64},
            },
        }()
        data_dict["x"] = collect(x)
        data_dict["y"] = collect(y)
        data_dict["func"] = convert(Array{Float64}, data_matrix)

        data_dict["xmin"] = xmin
        data_dict["xmax"] = xmax
        data_dict["ymin"] = ymin
        data_dict["ymax"] = ymax

        save(filepath, data_dict)
        @info "$(filename) created and exported!"
    end

    data_matrix[data_matrix.<1e-299] .= 1e-300

    function scaler(x::Real, scale::Symbol)
        if scale == :linear
            return x
        elseif scale == :log
            return log10(x)
        else
            error("Scale <$scale> is not supported")
        end
    end

    function un_scaler(x::Real, scale::Symbol)
        if scale == :linear
            return x
        elseif scale == :log
            return 10 .^ x
        else
            error("Scale <$scale> is not supported")
        end
    end


    knots = (scaler.(x, scale_x), scaler.(y, scale_y))
    f_matrix = scaler.(data_matrix, scale_f)

    itp = LinearInterpolation(knots, f_matrix , extrapolation_bc = extrapolation_bc())

    function func_interp(x::Real, y::Real)
        return un_scaler(itp(scaler(x, scale_x), scaler(y, scale_y)), scale_f)::Float64
    end

    return func_interp::Function
end


function create_interpolation_3D(
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
    scale_x=:linear,
    scale_y=:linear,
    scale_z=:linear,
    scale_f=:linear,
    extrapolation_bc = Throw,
    kwargs...,
)

    func_name = nameof(func)

    arg(x,y,z) = func(x, y, z, args...; kwargs...)

    if isnothing(jld_base_path)  
        base_path = pwd()
    else
        base_path = jld_base_path
    end

    if !(ispath(base_path))
        mkdir(base_path)
    end

    if isnothing(custom_name)
        filename = "$(func_name)_data.jld"
    else
        filename = "$(custom_name)_data.jld"
    end
    
    filepath = "$(base_path)/$(filename)"

    function load_file()
        if isfile(filepath)
            data = load(filepath)
            if data["xmin"] == xmin && data["xmax"] == xmax && data["ymin"] == ymin && data["ymax"] == ymax && data["zmin"] == zmin && data["zmax"] == zmax 
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
        if scale_x == :linear
            x = range(xmin, xmax, length = npoints_x)
        elseif scale_x == :log
            x = 10 .^ range(log10(xmin), log10(xmax), length = npoints_x)
        else
            error("X scale $scale_x not supported")
        end

        if scale_y == :linear
            y = range(ymin, ymax, length = npoints_y)
        elseif scale_y == :log
            y = 10 .^ range(log10(ymin), log10(ymax), length = npoints_y)
        else
            error("Y scale $scale_y not supported")
        end

        if scale_z == :linear
            z = range(zmin, zmax, length = npoints_z)
        elseif scale_z == :log
            z = 10 .^ range(log10(zmin), log10(zmax), length = npoints_z)
        else
            error("Z scale $scale_z not supported")
        end

        if isnothing(custom_name)
            @info "Computing $(func_name) Interpolation"
        else
            @info "Computing $(custom_name) Interpolation"
        end
        
        p = Progress(Int(npoints_z*npoints_y))
        update!(p, 0)

        data_matrix = zeros(npoints_x, npoints_y, npoints_z)

        @threads for k = 1:npoints_z
            for j = 1:npoints_y
                for i = 1:npoints_x
                    data_matrix[i] = arg(x[i], y[j], z[k])
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
            },
        }()
        data_dict["x"] = collect(x)
        data_dict["y"] = collect(y)
        data_dict["z"] = collect(z)
        data_dict["func"] = convert(Array{Float64}, data_matrix)

        data_dict["xmin"] = xmin
        data_dict["xmax"] = xmax
        data_dict["ymin"] = ymin
        data_dict["ymax"] = ymax
        data_dict["zmin"] = zmin
        data_dict["zmax"] = zmax

        save(filepath, data_dict)
        @info "$(filename) created and exported!"
    end

    data_matrix[data_matrix.<1e-299] .= 1e-300

    function scaler(x::Real, scale::Symbol)
        if scale == :linear
            return x
        elseif scale == :log
            return log10(x)
        else
            error("Scale <$scale> is not supported")
        end
    end

    function un_scaler(x::Real, scale::Symbol)
        if scale == :linear
            return x
        elseif scale == :log
            return 10 .^ x
        else
            error("Scale <$scale> is not supported")
        end
    end


    knots = (scaler.(x, scale_x), scaler.(y, scale_y), scaler.(z, scale_z))
    f_matrix = scaler.(data_matrix, scale_f)

    itp = LinearInterpolation(knots, f_matrix , extrapolation_bc = extrapolation_bc())

    function func_interp(x::Real, y::Real)
        return un_scaler(itp(scaler(x, scale_x), scaler(y, scale_y), scaler.(z, scale_z)), scale_f)::Float64
    end

    return func_interp::Function
end

#TODO
# function create_interpolation_ND(
#     func::Function;
#     jld_base_path = nothing,
#     custom_name = nothing,
#     xmin::Vector,
#     xmax::Vector,
#     npoints::Vector{Int},
#     scale_x::Union{Symbol, Vector{Symbol}} = :linear,
#     scale_f = :linear,
#     extrapolation_bc = Throw,
#     kwargs...,
# )

#     func_name = nameof(func)

#     arg(x) = func(x...; kwargs...)

#     if isnothing(jld_base_path)  
#         base_path = pwd()
#     else
#         base_path = jld_base_path
#     end

#     if !(ispath(base_path))
#         mkdir(base_path)
#     end

#     if isnothing(custom_name)
#         filename = "$(func_name)_data.jld"
#     else
#         filename = "$(custom_name)_data.jld"
#     end
    
#     filepath = "$(base_path)/$(filename)"

#     function load_file()
#         if isfile(filepath)
#             data = load(filepath)
#             if data["xmin"] == xmin && data["xmax"] == xmax 
#                 return true
#             else
#                 return false
#             end
#         else
#             return false
#         end
#     end

#     if load_file()
#         data = load(filepath)
#         @info "$(filename) loaded!"
#         x = data["x"]

#         data_matrix = data["func"]
#         #TODO ->
#     else
#         if scale_x == :linear
#             x = range(xmin, xmax, length = npoints)
#         elseif scale_x == :log
#             x = 10 .^ range(log10(xmin), log10(xmax), length = npoints)
#         else
#             error("X scale $scale_x not supported")
#         end

#         if isnothing(custom_name)
#             @info "Computing $(func_name) Interpolation"
#         else
#             @info "Computing $(custom_name) Interpolation"
#         end
        
#         p = Progress(Int(npoints))
#         update!(p, 0)

#         data_matrix = zeros(npoints)

#         @threads for i = 1:npoints

#             data_matrix[i] = arg(x[i])

#             next!(p)
#         end

#         data_dict = Dict{
#             String,
#             Union{
#                 String,
#                 Real,
#                 Array{Float64},
#             },
#         }()
#         data_dict["x"] = collect(x)
#         data_dict["func"] = convert(Array{Float64}, data_matrix)

#         data_dict["xmin"] = xmin
#         data_dict["xmax"] = xmax

#         save(filepath, data_dict)
#         @info "$(filename) created and exported!"
#     end

#     data_matrix[data_matrix.<1e-299] .= 1e-300

#     function scaler(x::Real, scale::Symbol)
#         if scale == :linear
#             return x
#         elseif scale == :log
#             return log10(x)
#         else
#             error("Scale <$scale> is not supported")
#         end
#     end

#     function un_scaler(x::Real, scale::Symbol)
#         if scale == :linear
#             return x
#         elseif scale == :log
#             return 10 .^ x
#         else
#             error("Scale <$scale> is not supported")
#         end
#     end


#     knots = (scaler.(x, scale_x), )
#     f_matrix = scaler.(data_matrix, scale_f)

#     itp = LinearInterpolation(knots, f_matrix , extrapolation_bc = extrapolation_bc())

#     function func_interp(x::Real)
#         return un_scaler(itp(scaler(x, scale_x)), scale_f)::Float64
#     end

#     return func_interp::Function
# end

end
