include("node.jl")

function tile_to_shape_along_axis(arr, target_shape, target_axis)
    # println("arr")
    # println(arr)
    # println("target_shape")
    # println(target_shape)
    # println("target_axis")
    # println(target_axis)
    if length(size(arr)) == 0
        return(repeat([arr], outer=target_shape))
    elseif length(size(arr)) == 1 && size(arr)[1] == target_shape[target_axis]
        repeating_shape = [i for i in target_shape]
        repeating_shape[target_axis] = 1
        arr_shape = ones(Int64, length(target_shape))
        arr_shape[target_axis] = length(arr)
        # Think I need to reshape arr to eb a matrix with 1 in all directions except the target exis
        return repeat(reshape(arr, Tuple(arr_shape)), outer=repeating_shape)
    else
        println("Not implemented")
    end
end

function tile_to_other_dist_along_axis_name(tiling_labeled_array::LabelledArray, target_array::LabelledArray)
    # println("tiling_labeled_array")
    # println(tiling_labeled_array)
    # println("target_array")
    # println(target_array)
    target_axis_label = tiling_labeled_array.axes_labels[1]
    target_axis_index = findfirst(==(target_axis_label), target_array.axes_labels)
    return LabelledArray(
        tile_to_shape_along_axis(tiling_labeled_array.array, size(target_array.array), target_axis_index),
        target_array.axes_labels
    )
end

# p_h1 = LabelledArray([0.6,0.3,0.1], ["h1"])

# p_v1_given_h1 = LabelledArray(
#     [0.4 0.8 0.9
#      0.6 0.2 0.1],
#      ["v1", "h1"]
# )

function other_axes_from_labeled_axes(labelled_array::LabelledArray, axis_label::String)
    return Tuple([i for i in 1:length(labelled_array.axes_labels) if labelled_array.axes_labels[i] != axis_label])
end

function variable_to_factor_messages(variable::Variable)
    # This needs to update the messsages in the factors from this variable
    for i in 1:length(variable.neighbours)
        # The i message needs to be excluded from the calculation
        new_message = ones(size(variable.incoming_messages[variable.neighbours[i].name]))
        for (key, incoming_message) in variable.incoming_messages
            if key != variable.neighbours[i].name
                new_message = new_message .* incoming_message
            end
        end
        variable.neighbours[i].incoming_messages[variable.name] = new_message
    end
end

function factor_to_variable_messages(factor::Factor)
    # This needs to update all the incoming messages of the connected variables
    for i in 1:length(factor.neighbours)
        factor_dist = copy(factor.data.array)
        neighbour_variable_names = [var.name for var in factor.neighbours if var.name != factor.neighbours[i].name]
        incoming_messages = [factor.incoming_messages[neighbour_name] for neighbour_name in neighbour_variable_names]
        # println(neighbour_variable_names)
        # println(length(incoming_messages))
        for j in 1:length(incoming_messages)
            # println(j)
            # println(incoming_messages[j])
            tiled_result = tile_to_other_dist_along_axis_name(LabelledArray(incoming_messages[j], [neighbour_variable_names[j]]), factor.data).array
            factor_dist = factor_dist .* tiled_result
        end
        other_axes = other_axes_from_labeled_axes(factor.data, factor.neighbours[i].name)
        value_to_squeeze = factor_dist
        for axis in other_axes
            value_to_squeeze = sum(value_to_squeeze, dims=axis)
        end
        factor.neighbours[i].incoming_messages[factor.name] = dropdims(value_to_squeeze; dims=Tuple(findall(size(value_to_squeeze) .== 1)))
    end
end

function marginal(variable::Variable)
    unnorm_p = :nothing
    for (key, val) in variable.incoming_messages
        if unnorm_p == :nothing
            unnorm_p = val
        else
            unnorm_p = unnorm_p .* val
        end
    end
    return unnorm_p ./ sum(unnorm_p)
end

function add_edge_between(variable::Variable, factor::Factor)
    push!(variable.neighbours, factor)
    push!(factor.neighbours, variable)
    factor.incoming_messages[variable.name] = 1.
    variable.incoming_messages[factor.name] = 1.
end