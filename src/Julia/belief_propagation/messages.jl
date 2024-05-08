using NaNMath, Base.Threads, Hadamard, FFTW
include("node.jl")

calculate_entropy(prob_dist) = -NaNMath.sum(prob_dist .* log2.(prob_dist))

normalisation_constant::Float64 = 1e6
addition_away_from_zero::Float64 = 1e-6

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

function variable_to_factor_messages(variable::Variable{AbsFactor}, damping_factor::Float64=1.0)
    # This needs to update the messsages in the factors from this variable
    # We need to make sure that the dist is the last factor added onto each variable
    # so that all values above them are discarded allowing us to have lots of different
    # for a variable without needing to specify that we don't need to send data to each one
    neighbours_to_include = ones(Bool, size(variable.incoming_messages)[1])
    for (i, neighbour) in enumerate(variable.neighbours)
        # println("Help")
        # println(i)
        # println(variable.neighbour_index_to_avoid)
        if variable.neighbour_index_to_avoid <= 0 || i < variable.neighbour_index_to_avoid
            # The i message needs to be excluded from the calculation
            neighbours_to_include[i] = false
            new_message = prod(variable.incoming_messages[neighbours_to_include, :], dims=1)[1, :]
            # Here normalise the output to be a prob dist.
            new_message ./= sum(new_message)
            new_message .*= normalisation_constant
            # println("Setting output: ", i, " ", new_message)
            neighbour.incoming_messages[variable.index_in_neighbours_neighbour[i]] = (damping_factor * new_message) .+ ((1 - damping_factor) * neighbour.incoming_messages[variable.index_in_neighbours_neighbour[i]])
            # neighbour.incoming_messages[variable.index_in_neighbours_neighbour[i]] ./= sum(neighbour.incoming_messages[variable.index_in_neighbours_neighbour[i]])
            neighbours_to_include[i] = true
        end
    end
end

function update_with_damping(variable::AbsVariable, damping_factor::Float64, message::Vector{Float64}, index_in_neighbours::Int64)
    if length(variable.incoming_messages[index_in_neighbours]) != length(message)
        variable.incoming_messages[index_in_neighbours, :] = message
    else
        variable.incoming_messages[index_in_neighbours, :] = (damping_factor * message) .+ ((1 - damping_factor) * variable.incoming_messages[index_in_neighbours])
        # variable.incoming_messages[index_in_neighbours, :] ./= sum(variable.incoming_messages[index_in_neighbours, :])
    end
end

function factor_to_variable_messages(factor::Factor{AbsVariable}, damping_factor::Float64=1.0)
    # This needs to update all the incoming messages of the connected variables
    tiled_incoming_messages = [tile_to_other_dist_along_axis_name(LabelledArray(factor.incoming_messages[i], [neighbour.name]), factor.data).array for (i, neighbour) in enumerate(factor.neighbours)]
    # println("None XOR")
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
        # message_out ./= sum(message_out)
        # neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i], :] = message_out
        update_with_damping(neighbour, damping_factor, message_out, factor.index_in_neighbours_neighbour[i])
    end
end

function factor_to_variable_messages(factor::XorFactor{AbsVariable}, damping_factor::Float64=1.0)
    # This needs to update all the incoming messages of the connected variables
    # This can make use of th fast walsh hadamard transform to be able to more efficiently compute the xor
    # println("XOR running")
    expected_size = 1 << factor.neighbours[1].number_of_bits
    transformed_incoming_messages = ones(length(factor.incoming_messages), expected_size)
    for i in eachindex(factor.incoming_messages)
        if length(factor.incoming_messages[i]) == expected_size
            transformed_incoming_messages[i, :] = fwht(factor.incoming_messages[i])
        else
            transformed_incoming_messages[i, :] = fwht(transformed_incoming_messages[i, :])
        end
    end
    neighbours_to_include = ones(Bool, size(factor.incoming_messages)[1])
    for (i, neighbour) in enumerate(factor.neighbours)
        neighbours_to_include[i] = false
        # The max with 0 is to ensure that we do not end up with any negative values beign passed around which can cause
        # issues in other places the negative value are all incredibly small when they start out anyway
        message_out = max.(0., ifwht(prod(transformed_incoming_messages[neighbours_to_include, :], dims=1)[1, :]))
        # if sum(message_out) < 0.125
        #     println("Possible problem")
        # end
        message_out ./= sum(message_out)
        message_out .*= normalisation_constant
        message_out .+= addition_away_from_zero
        # neighbour.incoming_messages[factor.index_in_neighbours_neighbour[i], :] = message_out
        update_with_damping(neighbour, damping_factor, message_out, factor.index_in_neighbours_neighbour[i])
        neighbours_to_include[i] = true
    end
end

function factor_to_variable_messages(factor::AddFactor{AbsVariable}, damping_factor::Float64=1.0)
    # This relies on making the updates to the incoming messages based on what we know about them
    # They need to be ordered as carry_in, a, b, output and then we can make use of a series of 
    # fourier transforms and theorems from that to produce the correct set of results for a particular
    # way of doing it

    # First enlarge all the messages to be the same size
    # println(factor.neighbours[4].number_of_bits)
    size_of_variables = 1 << (factor.neighbours[4].number_of_bits)
    size_of_incoming_variables = size_of_variables ÷ 2
    carry_in_message = zeros(size_of_variables)
    if length(factor.incoming_messages[1]) == 1
        carry_in_message[1:2] .= 1.0
    else
        # println(factor.incoming_messages[1])
        carry_in_message[1:2] = factor.incoming_messages[1]
    end

    a_message = zeros(size_of_variables)
    if length(factor.incoming_messages[2]) == 1
        a_message[1:size_of_incoming_variables] .= 1.0
    else
        a_message[1:size_of_incoming_variables] = factor.incoming_messages[2]
    end

    b_message = zeros(size_of_variables)
    if length(factor.incoming_messages[3]) == 1
        b_message[1:size_of_incoming_variables] .= 1.0
    else
        b_message[1:size_of_incoming_variables] = factor.incoming_messages[3]
    end

    out_message = zeros(size_of_variables)
    if length(factor.incoming_messages[4]) == 1
        out_message[:] .= 1.0
    else
        out_message = factor.incoming_messages[4]
    end

    # println(carry_in_message)
    # println(a_message)
    # println(b_message)
    # println(out_message)

    C_IN = factor.fft_plan * carry_in_message
    A = factor.fft_plan * a_message
    B = factor.fft_plan * b_message
    OUT = factor.fft_plan * out_message

    t_c_in = real(factor.ifft_plan * (conj.(A .* B) .* OUT))
    # t_c_in = real(ifft(conj.(A) .* conj(B) .* OUT))
    t_c_in = max.(0.0, t_c_in[1:2])
    t_c_in ./= sum(t_c_in)
    t_c_in .*= normalisation_constant
    t_c_in .+= addition_away_from_zero
    # println(t_c_in)


    t_a = real(factor.ifft_plan * (conj.(C_IN .* B) .* OUT))
    # t_a = real(ifft(conj.(C_IN) .* conj(B) .* OUT))
    t_a = max.(0.0, t_a[1:size_of_incoming_variables])
    t_a ./= sum(t_a)
    t_a .*= normalisation_constant
    t_a .+= addition_away_from_zero
    # println(t_a)

    t_b = real(factor.ifft_plan * (conj.(A .* C_IN) .* OUT))
    # t_b = real(ifft(conj.(A) .* conj(C_IN) .* OUT))
    t_b = max.(0.0, t_b[1:size_of_incoming_variables])
    t_b ./= sum(t_b)
    t_b .*= normalisation_constant
    t_b .+= addition_away_from_zero
    # println(t_b)


    t_out = max.(0.0, real(factor.ifft_plan * (A .* B .* C_IN)))
    t_out ./= sum(t_out)
    t_out .*= normalisation_constant
    t_out .+= addition_away_from_zero
    # println(t_out)

    # Need to update the messages going out of from this factor to what has just come so shrink them and normalise them
    update_with_damping(factor.neighbours[1], damping_factor, t_c_in, factor.index_in_neighbours_neighbour[1])
    update_with_damping(factor.neighbours[2], damping_factor, t_a, factor.index_in_neighbours_neighbour[2])
    update_with_damping(factor.neighbours[3], damping_factor, t_b, factor.index_in_neighbours_neighbour[3])
    update_with_damping(factor.neighbours[4], damping_factor, t_out, factor.index_in_neighbours_neighbour[4])
end

function factor_to_variable_messages(factor::MarginaliseTopBitsFactor{AbsVariable}, damping_factor::Float64=1.0)
    # The first neighbour is the larger version and the second neighbour the smaller
    # println("Top bits")
    number_of_bits_large_variable = factor.neighbours[1].number_of_bits
    number_of_bits_smaller_variable = factor.neighbours[2].number_of_bits
    size_of_small_in_large = 1 << (number_of_bits_large_variable - number_of_bits_smaller_variable)
    # println(number_of_bits_large_variable)
    # println(number_of_bits_smaller_variable)
    # println(size_of_small_in_large)

    # We only need to consider the other neighbours incoming message for the out going message
    if length(factor.incoming_messages[1]) == 1
        small_outgoing_message = ones(1 << number_of_bits_smaller_variable) ./ (1 << number_of_bits_smaller_variable)
    else
        small_outgoing_message = collect(Iterators.map(sum, Iterators.partition(factor.incoming_messages[1], size_of_small_in_large)))
    end
    if length(factor.incoming_messages[2]) == 1
        large_outgoing_message = ones(1 << number_of_bits_large_variable)
    else
        large_outgoing_message = repeat(factor.incoming_messages[2], inner=size_of_small_in_large)
    end
    # large_outgoing_message ./= sum(large_outgoing_message)

    # Need to update the messages going out of from this factor to what has just come so shrink them and normalise them
    update_with_damping(factor.neighbours[1], damping_factor, large_outgoing_message, factor.index_in_neighbours_neighbour[1])
    update_with_damping(factor.neighbours[2], damping_factor, small_outgoing_message, factor.index_in_neighbours_neighbour[2])
end

function marginalise_bottom(values::Vector{Float64}, number_of_values::Int64)
    output = zeros(number_of_values)
    for i in 1:number_of_values
        output[i] = sum(values[i:number_of_values:end])
    end
    return output
end

function factor_to_variable_messages(factor::MarginaliseBottomBitsFactor{AbsVariable}, damping_factor::Float64=1.0)
    # The first neighbour is the larger version and the second neighbour the smaller
    # println("Bottom bits")
    number_of_bits_large_variable = factor.neighbours[1].number_of_bits
    number_of_bits_smaller_variable = factor.neighbours[2].number_of_bits
    size_of_small_in_large = 1 << (number_of_bits_large_variable - number_of_bits_smaller_variable)
    # println(number_of_bits_large_variable)
    # println(number_of_bits_smaller_variable)
    # println(size_of_small_in_large)
    # We only need to consider the other neighbours incoming message for the out going message
    if length(factor.incoming_messages[1]) == 1
        small_outgoing_message = ones(1 << number_of_bits_smaller_variable) ./ (1 << number_of_bits_smaller_variable)
    else
        small_outgoing_message = marginalise_bottom(factor.incoming_messages[1], 1 << number_of_bits_smaller_variable)
    end
    if length(factor.incoming_messages[2]) == 1
        large_outgoing_message = ones(1 << number_of_bits_large_variable)
    else
        large_outgoing_message = repeat(factor.incoming_messages[2], outer=size_of_small_in_large)
    end
    # large_outgoing_message ./= sum(large_outgoing_message)

    # Need to update the messages going out of from this factor to what has just come so shrink them and normalise them
    update_with_damping(factor.neighbours[1], damping_factor, large_outgoing_message, factor.index_in_neighbours_neighbour[1])
    update_with_damping(factor.neighbours[2], damping_factor, small_outgoing_message, factor.index_in_neighbours_neighbour[2])
end

function factor_to_variable_messages(factor::RotateFactor{AbsVariable}, damping_factor::Float64=1.0)
    # The first neighbour is the left, second neighbour is the right, and third is output
    # println("Rotate")
    number_of_bits_cluster = factor.neighbours[1].number_of_bits
    shift_amount = factor.bits_to_rotate_by
    # println(number_of_bits_cluster)
    # println(shift_amount)

    # We only need to consider the other neighbours incoming message for the out going message
    if length(factor.incoming_messages[1]) == 1
        left_incoming_message = ones(1 << number_of_bits_cluster)
    else
        left_incoming_message = repeat(marginalise_bottom(factor.incoming_messages[1], 1 << (number_of_bits_cluster - shift_amount)), inner=1 << shift_amount)
    end
    # println(left_incoming_message)
    if length(factor.incoming_messages[2]) == 1
        right_incoming_message = ones(1 << number_of_bits_cluster)
    else
        right_incoming_message = repeat(collect(Iterators.map(sum, Iterators.partition(factor.incoming_messages[2], 1 << (number_of_bits_cluster - shift_amount)))), outer=1 << (number_of_bits_cluster - shift_amount))
    end
    # println(right_incoming_message)
    if length(factor.incoming_messages[3]) == 1
        output_incoming_message = ones(1 << number_of_bits_cluster)
    else
        output_incoming_message = factor.incoming_messages[3]
    end
    # println(output_incoming_message)
    # Think I need to marginalise out the top bits and convert to bottom bits
    # println(collect(Iterators.map(sum, Iterators.partition(right_incoming_message .* output_incoming_message, 1 << (number_of_bits_cluster - shift_amount)))))
    # println(repeat(collect(Iterators.map(sum, Iterators.partition(right_incoming_message .* output_incoming_message, 1 << (number_of_bits_cluster - shift_amount)))), outer=1 << shift_amount))
    left_outgoing_message = repeat(collect(Iterators.map(sum, Iterators.partition(right_incoming_message .* output_incoming_message, 1 << shift_amount))), outer=1 << shift_amount)
    # left_outgoing_message ./= sum(left_outgoing_message)

    # Similar need to marginalise the bottom bits and covert to top bits
    right_outgoing_message = repeat(marginalise_bottom(left_incoming_message .* output_incoming_message, 1 << shift_amount), inner=1 << (number_of_bits_cluster - shift_amount))
    # right_outgoing_message ./= sum(right_outgoing_message)

    output_outgoing_message = left_incoming_message .* right_incoming_message
    # output_outgoing_message ./= sum(output_outgoing_message)

    # Need to update the messages going out of from this factor to what has just come so shrink them and normalise them
    update_with_damping(factor.neighbours[1], damping_factor, left_outgoing_message, factor.index_in_neighbours_neighbour[1])
    update_with_damping(factor.neighbours[2], damping_factor, right_outgoing_message, factor.index_in_neighbours_neighbour[2])
    update_with_damping(factor.neighbours[3], damping_factor, output_outgoing_message, factor.index_in_neighbours_neighbour[3])
end

function marginal(variable::AbsVariable)
    unnorm_p = prod(variable.incoming_messages, dims=1)[1, :]
    total_unorm_p = NaNMath.sum(unnorm_p)
    if isnan(total_unorm_p)
        return zeros(length(unnorm_p))
    end
    return total_unorm_p > 0 ? unnorm_p ./ total_unorm_p : unnorm_p
end

function add_edge_between(variable::AbsVariable, factor::AbsFactor)
    push!(variable.neighbours, factor)
    push!(factor.neighbours, variable)

    push!(factor.incoming_messages, [1.0])
    variable.incoming_messages = ones(size(variable.incoming_messages)[1] + 1, 1 << variable.number_of_bits) ./ 1 << variable.number_of_bits
    # push!(variable.incoming_messages, [1.])

    push!(factor.index_in_neighbours_neighbour, length(variable.neighbours))
    push!(variable.index_in_neighbours_neighbour, length(factor.neighbours))
end


function set_variable_to_value(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    variable_name_with_version::String,
    value::UInt32,
    number_of_bits_per_cluster::Int64,
    run_number::Int64
)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    for i in 1:number_of_clusters
        cur_var_name = string(variable_name_with_version, "_", i, "_", run_number)
        cur_dist_name = string("f_", cur_var_name, "_dist")
        dist_table = zeros(1 << number_of_bits_per_cluster)
        # Calculate what value these bits should have
        cur_cluster_value = (value & (((1 << number_of_bits_per_cluster) - 1)) << (number_of_bits_per_cluster * (i - 1))) >> (number_of_bits_per_cluster * (i - 1))
        dist_table[cur_cluster_value+1] = 1.0
        factors[cur_dist_name] = Factor{AbsVariable}(cur_dist_name, LabelledArray(dist_table, [cur_var_name]))
        add_edge_between(variables[cur_var_name], factors[cur_dist_name])
        variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
    end
end

function read_most_likely_value_from_variable(variables::Dict{String,AbsVariable},
    variable_name_with_version::String,
    number_of_bits_per_cluster::Int64,
    run_number::Int64)

    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    value::UInt32 = 0
    for i in 1:number_of_clusters
        cur_var_name = string(variable_name_with_version, "_", i, "_", run_number)
        dist = marginal(variables[cur_var_name])
        # println(cur_var_name,": ", dist)
        value += (findmax(dist)[2] - 1) << (number_of_bits_per_cluster * (i - 1))
    end
    return value
end

function update_all_entropies(variables::Dict{String,AbsVariable},
    all_var_names::Vector{String})
    Threads.@threads for var_name in all_var_names
        update_variable_entropy(variables[var_name])
    end
end

function update_variable_entropy(variable::AbsVariable)
    variable.previous_entropy = variable.current_entropy
    variable.current_entropy = calculate_entropy(marginal(variable))
end

function total_entropy_of_graph(variables::Dict{String,AbsVariable})
    tot_ent = 0.0
    for (i, j) in variables
        tot_ent += j.current_entropy
    end
    return tot_ent
end

function belief_propagate_forwards_and_back_through_graph(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    variables_by_round::Vector{Set{String}},
    factors_by_round::Vector{Set{String}},
    times_per_round::Int64,
    damping_factor::Float64)
    # GO forward first
    for (i, vars_for_round) in enumerate(variables_by_round)
        for j in 1:times_per_round
            for var_name in vars_for_round
                variable_to_factor_messages(variables[var_name], damping_factor)
            end
            for fact_name in factors_by_round[i]
                factor_to_variable_messages(factors[fact_name], damping_factor)
            end
        end
    end

    # # Search for a variable which has an invalid distribution to see what we end up with
    # for var in keys(variables)
    #     if abs(sum(marginal(variables[var])) - 1) >= 1e-2
    #         println(var)
    #     end
    # end

    # Then go backwards
    for i in length(variables_by_round):-1:1
        for j in 1:times_per_round
            for var_name in variables_by_round[i]
                variable_to_factor_messages(variables[var_name], damping_factor)
            end
            for fact_name in factors_by_round[i]
                factor_to_variable_messages(factors[fact_name], damping_factor)
            end
        end
    end

    # for var in keys(variables)
    #     if abs(sum(marginal(variables[var])) - 1) >= 1e-2
    #         println(var)
    #     end
    # end
end

# From looking at a few of the animations of how entropy is changing in the graph over time
# it becomes clear that the most changes happen at the ends of the execution therefore 
# I think a good apporach can be to do a set number of full iterations through the graph (like 10)
# before just doing the first and final x rounds of the propagation to hopefully improve performances and 
# then might be able to just return values from that once they are no longer changing or after a certain number of rounds
# for that part