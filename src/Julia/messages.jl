include("node.jl")

messages = Dict()

function _tile_to_shape_along_axis(arr, target_shape, target_axis)

end

function _tile_to_other_dist_along_axis_name(tiling_labeled_array::LabelledArray, target_array::LabelledArray)
    target_axis_label = tiling_labeled_array.axes_labels[1]
    return LabelledArray(
        _tile_to_shape_along_axis(tiling_labeled_array.array, size(target_array.array), target_axis_label), #The final part really needs to be changed to give the index rather than the label
        target_array.axes_names
    )
end

function _variable_to_factor_messages(variable::Variable, factor::Factor)
    incoming_messages = [] # This needs to be the messages excluding those involving factor
    return prod(incoming_messages, axes=1)
end

function _factor_to_variable_messages(factor::Factor, variable::Variable)
    
end