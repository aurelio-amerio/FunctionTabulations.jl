function _load_file_2D(filepath, xmin, xmax, ymin, ymax, npoints_x, npoints_y, metadata, metadata_validation_fn)
    if isfile(filepath)
        data = load(filepath)
        if data["xmin"] == xmin && data["xmax"] == xmax && data["ymin"] == ymin && data["ymax"] == ymax && data["npoints_x"] == npoints_x && data["npoints_y"] == npoints_y
            if metadata === nothing
                return true
            else
                return check_metadata(metadata, data["metadata"], metadata_validation_fn)
            end
        else
            @warn "The provided tabulation range or npoints does not match the data stored in the tabulation file. The tabulation will be recomputed."

            return false
        end
    else
        return false
    end
end

function _compute_matrix_and_save_2D(func, arg, x, y, x_scale, y_scale, xmin, xmax, ymin, ymax, npoints_x, npoints_y, metadata, filepath, filename)
    p = Progress(Int(npoints_y))
    update!(p, 0)

    data_matrix = zeros(npoints_x, npoints_y)

    for j = 1:npoints_y
        @threads for i = 1:npoints_x
            data_matrix[i, j] = arg(un_scaler(x[i], x_scale), un_scaler(y[j], y_scale))
        end
        next!(p)
    end

    data_dict = Dict{String,Any}()
    sha = compute_SHA(func)
    data_dict["x"] = x
    data_dict["y"] = y
    data_dict["func"] = convert(Array{Float64}, data_matrix)
    data_dict["sha"] = sha
    data_dict["xmin"] = xmin
    data_dict["xmax"] = xmax
    data_dict["ymin"] = ymin
    data_dict["ymax"] = ymax
    data_dict["npoints_x"] = npoints_x
    data_dict["npoints_y"] = npoints_y
    data_dict["metadata"] = metadata

    save(filepath, data_dict)
    @info "$(filename) created and exported!"
    return data_matrix
end

function _get_interp_fn_2D(x, y, data_matrix, x_scale, y_scale, f_scale, x_units, y_units, f_units, interpolation_type, extrapolation_bc)
    if f_scale == :log10
        data_matrix[data_matrix.<=1e-300] .= 1e-300
    end

    knots = (x, y)
    f_matrix = scaler.(data_matrix, f_scale)

    if interpolation_type == :linear
        itp = linear_interpolation(knots, f_matrix, extrapolation_bc=extrapolation_bc())
    elseif interpolation_type == :cubic
        itp = cubic_spline_interpolation(knots, f_matrix, extrapolation_bc=extrapolation_bc())
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end

    function func_interp(x::Real, y::Real)
        return un_scaler(itp(scaler(x, x_scale), scaler(y, y_scale)), f_scale)::Float64
    end

    return wrap_function_2D_add_units(func_interp, x_units, y_units, f_units)
end