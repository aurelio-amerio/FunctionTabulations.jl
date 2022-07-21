function _load_file_1D(filepath, xmin, xmax, npoints, metadata, metadata_validation_fn)
    if isfile(filepath)
        data = load(filepath)
        if data["xmin"] == xmin && data["xmax"] == xmax && data["npoints"] == npoints
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

function _compute_matrix_and_save_1D(func, arg, x, x_scale, xmin, xmax, npoints, metadata, filepath, filename)
    p = Progress(Int(npoints))
    update!(p, 0)

    data_matrix = zeros(npoints)

    @threads for i = 1:npoints

        data_matrix[i] = arg(un_scaler(x[i], x_scale))

        next!(p)
    end

    data_dict = Dict{String,Any}()
    sha = compute_SHA(func)
    data_dict["x"] = x
    data_dict["func"] = convert(Array, data_matrix)
    data_dict["sha"] = sha
    data_dict["xmin"] = xmin
    data_dict["xmax"] = xmax
    data_dict["npoints"] = npoints
    data_dict["metadata"] = metadata

    save(filepath, data_dict)
    @info "$(filename) created and exported!"
    return data_matrix
end

function _get_interp_fn_1D(x, data_matrix, x_scale, f_scale, x_units, f_units, interpolation_type, extrapolation_bc)
    if f_scale == :log10
        data_matrix[data_matrix.<=1e-300] .= 1e-300
    end

    knots = (x,)
    f_matrix = scaler.(data_matrix, f_scale)

    if interpolation_type == :linear
        itp = LinearInterpolation(knots, f_matrix, extrapolation_bc=extrapolation_bc())
    elseif interpolation_type == :cubic
        itp = CubicSplineInterpolation(knots, f_matrix, extrapolation_bc=extrapolation_bc())
    else
        throw(ArgumentError("$interpolation_type is not a valid interpolation type"))
    end


    function func_interp(x::Real)
        return un_scaler(itp(scaler(x, x_scale)), f_scale)::Float64
    end

    return wrap_function_1D_add_units(func_interp, x_units, f_units)
end