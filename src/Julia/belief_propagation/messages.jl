using NaNMath, Base.Threads
include("node.jl")

calculate_entropy(prob_dist) = -NaNMath.sum(prob_dist .* log2.(prob_dist))

function tile_to_shape_along_axis(arr::Float64, target_shape::Tuple, target_axis::Int64)
    return fill(arr, target_shape)
end

function tile_to_shape_along_axis(arr::Vector{Float64}, target_shape::Tuple, target_axis::Int64)
    if length(arr) == 1
        return fill(arr[1], target_shape)
    elseif length(arr) == target_shape[target_axis]
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
    target_axis_label = tiling_labeled_array.axes_labels[1]
    target_axis_index = findfirst(==(target_axis_label), target_array.axes_labels)
    return LabelledArray(
        tile_to_shape_along_axis(tiling_labeled_array.array, size(target_array.array), target_axis_index),
        target_array.axes_labels
    )
end

function other_axes_from_labeled_axes(labelled_array::LabelledArray, axis_label::String)
    return [i for i in 1:length(labelled_array.axes_labels) if labelled_array.axes_labels[i] != axis_label]
end

damping_factor = 1.

function variable_to_factor_messages(variable::Variable{Factor})
    # This needs to update the messsages in the factors from this variable
    neighbours_to_include = ones(Bool, size(variable.incoming_messages)[1])
    for (i, neighbour) in enumerate(variable.neighbours)
        if i != variable.neighbour_index_to_avoid
            # The i message needs to be excluded from the calculation
            neighbours_to_include[i] = false
            new_message = prod(variable.incoming_messages[neighbours_to_include, :], dims=1)[1, :]
            # Here normalise the output to be a prob dist.
            new_message ./= sum(new_message)
            neighbour.incoming_messages[variable.index_in_neighbours_neighbour[i]] = new_message#damping_factor * new_message .+ (1 - damping_factor) * neighbour.incoming_messages[variable.index_in_neighbours_neighbour[i]]
            neighbours_to_include[i] = true
        end
    end
end

function factor_to_variable_messages(factor::Factor{Variable})
    # This needs to update all the incoming messages of the connected variables
    tiled_incoming_messages = [tile_to_other_dist_along_axis_name(LabelledArray(factor.incoming_messages[i], [neighbour.name]), factor.data).array for (i, neighbour) in enumerate(factor.neighbours)]
    for (i, neighbour) in enumerate(factor.neighbours)
        factor_dist = copy(factor.data.array)
        for (j, tiled_incoming_message) in enumerate(tiled_incoming_messages)
            if i != j
                factor_dist .*= tiled_incoming_messages[j]
            end
        end
        other_axes = other_axes_from_labeled_axes(factor.data, factor.neighbours[i].name)
        value_to_squeeze = sum(factor_dist; dims=other_axes)
        message_out = dropdims(value_to_squeeze; dims=Tuple(other_axes))
        message_out ./= sum(message_out)
        neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i], :] = message_out
        # if length(neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i]]) != length(message_out)
        #     neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i], :] = message_out
        # else
        #     neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i], :] = damping_factor * message_out .+ (1 - damping_factor) * neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i]]
        # end
    end
end

function marginal(variable::Variable{Factor})
    unnorm_p = prod(variable.incoming_messages, dims=1)[1,:]
    total_unorm_p = NaNMath.sum(unnorm_p)
    if isnan(total_unorm_p)
        return zeros(length(unnorm_p))
    end
    return total_unorm_p > 0 ? unnorm_p ./ total_unorm_p : unnorm_p
end

function add_edge_between(variable::Variable{Factor}, factor::Factor{Variable})
    push!(variable.neighbours, factor)
    push!(factor.neighbours, variable)
    
    push!(factor.incoming_messages, [1.])
    variable.incoming_messages = ones(size(variable.incoming_messages)[1] + 1, 1 << variable.number_of_bits) ./ 1 << variable.number_of_bits
    # push!(variable.incoming_messages, [1.])
    
    push!(factor.index_in_neighbours_neighbour, length(variable.neighbours))
    push!(variable.index_in_neighbours_neighbour, length(factor.neighbours))
end


function set_variable_to_value(variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
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
        factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(dist_table, [cur_var_name]))
        add_edge_between(variables[cur_var_name], factors[cur_dist_name])
        variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
    end
end

function read_most_likely_value_from_variable(variables::Dict{String, Variable{Factor}},
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

function update_all_entropies(variables::Dict{String, Variable{Factor}},
    all_var_names::Vector{String})
    Threads.@threads for var_name in all_var_names
        update_variable_entropy(variables[var_name])
    end
end

function update_variable_entropy(variable::Variable{Factor})
    variable.previous_entropy = variable.current_entropy
    variable.current_entropy = calculate_entropy(marginal(variable))
end

function total_entropy_of_graph(variables::Dict{String, Variable{Factor}})
    tot_ent = 0.
    for (i,j) in variables
        tot_ent += j.current_entropy
    end
    return tot_ent
end

function belief_propagate_forwards_and_back_through_graph(variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    variables_by_round::Vector{Set{String}},
    factors_by_round::Vector{Set{String}},
    times_per_round::Int64)
    # GO forward first
    for (i, vars_for_round) in enumerate(variables_by_round)
        for j in 1:times_per_round
            for var_name in vars_for_round
                variable_to_factor_messages(variables[var_name])
            end
            for fact_name in factors_by_round[i]
                factor_to_variable_messages(factors[fact_name])
            end
        end
    end
    # Then go backwards
    for i in length(variables_by_round):-1:1
        for j in 1:times_per_round
            for var_name in variables_by_round[i]
                variable_to_factor_messages(variables[var_name])
            end
            for fact_name in factors_by_round[i]
                factor_to_variable_messages(factors[fact_name])
            end
        end
    end
end