include("node.jl")

messages = Dict()

function tile_to_shape_along_axis(arr, target_shape, target_axis)
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


function variable_to_factor_messages(variable::Variable)
    # This needs to update the messsages in the factors from this variable
    for i in 1:length(variable.neighbours)
        # The i message needs to be excluded from the calculation
        new_message = prod(variable.incoming_messages[begin:end .!=i], axes=1)
        index_in_factor = findfirst(==(variable), variable.neighbours[i].neighbours)
        variable.neighbours[i].incoming_messages[index_in_factor] = new_message
    end
end

function factor_to_variable_messages(factor::Factor)
    # This needs to update all the incoming messages of the connected variables
    for i in 1:length(factor.neighbours)
        factor_dist = copy(factor.data.array)
        incoming_messages = factor.incoming_messages[begin:end .!=i]
        neighbour_variable_names = factor.neighbours[begin:end .!=i].name
        tiled_results = tile_to_other_dist_along_axis_name.(LabelledArray.(incoming_messages, neighbour_variable_names), factor.data).array
        # Need to have the product of the above
        for tiled_res in tiled_results
            factor_dist = factor_dist .* tiled_res
        end
        other_axes = other_axes_from_labeled_axes(factor.data, factor.neighbours[i].name)
        value_to_squeeze = sum(factor_dist, axes=other_axes)
        index_in_variable = findfirst(==(factor), factor.neighbours[i].neighbours)
        factor.neighbours[i].incoming_messages[index_in_variable] = dropdims(value_to_squeeze; dims=Tuple(findall(size(value_to_squeeze) .== 1)))
    end
end