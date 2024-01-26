using NaNMath
include("node.jl")

function tile_to_shape_along_axis(arr::Float64, target_shape::Tuple, target_axis::Int64)
    return fill(arr, target_shape)
end

function tile_to_shape_along_axis(arr::Vector{Float64}, target_shape::Tuple, target_axis::Int64)
    # println("arr")
    # println(arr)
    # println("target_shape")
    # println(target_shape)
    # println("target_axis")
    # println(target_axis)
    if length(arr) == target_shape[target_axis]
        repeating_shape = [i for i in target_shape]
        repeating_shape[target_axis] = 1
        arr_shape = ones(Int64, length(target_shape))
        arr_shape[target_axis] = length(arr)
        # Think I need to reshape arr to be a matrix with 1 in all directions except the target eaxis
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
        for j in 1:length(neighbour_variable_names)
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
        # Normalise the message, this should hopefully stop values spinning off to infinity but I am not entirely sure if it is valid
        # Haven't done this on the variable to factor messages but this should hopefully be enough to keep it first_number_cluster_shifts
        # it might also be better to have it work by dividing by a power of 2 near the value so that the division does not really need to be done
        # because that would allow it to stay still contained
        message_out = dropdims(value_to_squeeze; dims=Tuple(findall(size(value_to_squeeze) .== 1)))
        message_out_sum = sum(message_out)
        message_out = message_out_sum > 0 ? message_out ./ message_out_sum : message_out
        if message_out_sum <= 0
            println(factor.name)
        end
        factor.neighbours[i].incoming_messages[factor.name] = message_out
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
    total_unorm_p = sum(unnorm_p)
    return total_unorm_p > 0 ? unnorm_p ./ total_unorm_p : unnorm_p
end

function add_edge_between(variable::Variable, factor::Factor)
    push!(variable.neighbours, factor)
    push!(factor.neighbours, variable)
    factor.incoming_messages[variable.name] = 1.
    variable.incoming_messages[factor.name] = 1.
end


function set_variable_to_value(variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    variable_name_with_version::String,
    value::UInt32,
    number_of_bits_per_cluster::Int64
    )
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    for i in 1:number_of_clusters
        cur_var_name = string(variable_name_with_version, "_", i)
        cur_dist_name = string("f_", cur_var_name, "_dist")
        dist_table = zeros(1 << number_of_bits_per_cluster)
        # Calculate what value these bits should have
        cur_cluster_value = (value & (((1 << number_of_bits_per_cluster) - 1)) << (number_of_bits_per_cluster * (i - 1))) >> (number_of_bits_per_cluster * (i - 1))
        dist_table[cur_cluster_value + 1] = 1.
        factors[cur_dist_name] = Factor(cur_dist_name, LabelledArray(dist_table, [cur_var_name]))
        add_edge_between(variables[cur_var_name], factors[cur_dist_name])
    end
end

function read_most_likely_value_from_variable(variables::Dict{String, Variable},
    variable_name_with_version::String,
    number_of_bits_per_cluster::Int64)

    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    value::UInt32 = 0
    for i in 1:number_of_clusters
        cur_var_name = string(variable_name_with_version, "_", i)
        dist = marginal(variables[cur_var_name])
        # println(cur_var_name,": ", dist)
        value += (findmax(dist)[2] - 1) << (number_of_bits_per_cluster * (i - 1))
    end
    return value
end

function total_entropy_of_graph(variables::Dict{String, Variable})
    tot_ent = 0.
    for (i,j) in variables
        prob_dist = marginal(j)
        tot_ent -= NaNMath.sum(prob_dist .* log2.(prob_dist))
    end
    return tot_ent
end
