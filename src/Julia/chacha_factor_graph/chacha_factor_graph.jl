include("../belief_propagation/node.jl")
include("../belief_propagation/messages.jl")

# This creates an entire factor graph for the ChaCha algorithm

# It works by having a numbering for the different locations in ChaCha algorithm as 
# 1  2  3  4 
# 5  6  7  8 
# 9  10 11 12
# 13 14 15 16

# Then going forward the variables will be labeled as {location}_{number_of_times_used}_{cluster_number}_{run_number},
# with the additional variables introduced by variables but not native in the algorithm (e.g. carry bits, temporary parts in adds)
# these will also be labeled by refering to the number of times they have been used alongside things such as location

function make_rotation_prob_table(number_of_bits_in_cluster::Int64, number_of_bits_to_shift::Int64)
    output = zeros(1 << number_of_bits_in_cluster, 1 << number_of_bits_in_cluster, 1 << number_of_bits_in_cluster)
    number_of_bits_from_right = number_of_bits_to_shift
    number_of_bits_from_left = number_of_bits_in_cluster - number_of_bits_from_right
    for input_left in 1:(1<<number_of_bits_in_cluster)
        for input_right in 1:(1<<number_of_bits_in_cluster)
            output[input_left, input_right, ((((input_left-1)&((1<<number_of_bits_from_left)-1))<<number_of_bits_from_right)|((((input_right-1)>>number_of_bits_from_left)&((1<<number_of_bits_from_right)-1))))+1] = 1.0
        end
    end
    return output
end

function make_xor_prob_table(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 1 << num_of_bits, 1 << num_of_bits)
    for input_a in 1:(1<<num_of_bits)
        for input_b in 1:(1<<num_of_bits)
            output[input_a, input_b, ((input_a-1)⊻(input_b-1))+1] = 1.0
        end
    end
    return output
end

function make_add_including_carry_prob_array(num_of_bits::Int64)
    output = zeros(2, 1 << num_of_bits, 1 << num_of_bits, 1 << (num_of_bits + 1))
    for carry_in in 1:2
        for input_a in 1:(1<<num_of_bits)
            for input_b in 1:(1<<num_of_bits)
                output[carry_in, input_a, input_b, input_a+input_b+carry_in-2] = 1.0
            end
        end
    end
    return output
end

function take_top_bit_prob_array(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 2)
    for val_in in 1:(1<<num_of_bits)
        output[val_in, (val_in>1<<(num_of_bits-1))+1] = 1.0
    end
    return output
end

function take_bottom_bits_prob_array(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 1 << (num_of_bits - 1))
    for val_in in 1:1<<num_of_bits
        output[val_in, ((val_in-1)%((1<<(num_of_bits-1))))+1] = 1.0
    end
    return output
end

function make_equality_prob_array(num_of_bits::Int64, number_of_inputs::Int64)
    output = zeros(Tuple([1 << num_of_bits for _ in 1:number_of_inputs]))
    for input_a in 1:(1<<num_of_bits)
        output[CartesianIndex(Tuple([input_a for _ in 1:number_of_inputs]))] = 1.0
    end
    return output
end

function bit_rotation_factor_graph!(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    input::Int64,
    bits_to_rotate_by::Int64,
    number_of_bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    number_of_operations::Dict{String,Int64},
    precalculated_prob_tables::Dict{String,Array{Float64}},
    round_variables::Set{String},
    round_factors::Set{String},
    run_number::Int64)
    if bits_to_rotate_by == 0
        return
    end
    # Need to add the logic for performing a left bitwise rotation of the value in the input
    # We know that the amount which is being shifted by is smaller than a cluster therefore we can effectively do a
    # left shift between adjacent clusters, e.g. the new value of cluster two takes the top bits_to_rotate_by of the 
    # lower cluster and (number_of_bits_per_cluster - bits_to_rotate_by) bottom bits of the higher cluster (e.g. further to the left)
    # this is handled inside of the probability table of the rotation currently

    rotation_cluster_prob_table = precalculated_prob_tables[string("rotation_cluster_", bits_to_rotate_by)]
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    for i in 1:number_of_clusters
        further_left_variable_name = string(input, "_", location_execution_counts[input], "_", i, "_", run_number)
        further_right_variable_name = string(input, "_", location_execution_counts[input], "_", i - 1 <= 0 ? i - 1 + number_of_clusters : i - 1, "_", run_number)
        output_variable_name = string(input, "_", location_execution_counts[input] + 1, "_", i, "_", run_number)
        factor_name = string("f_rot_", number_of_operations["rot"], "_", i, "_", run_number)
        push!(round_variables, further_left_variable_name)
        push!(round_variables, further_right_variable_name)
        push!(round_variables, output_variable_name)
        push!(round_factors, factor_name)

        variables[output_variable_name] = Variable{AbsFactor}(output_variable_name, number_of_bits_per_cluster)
        factors[factor_name] = Factor{AbsVariable}(factor_name, LabelledArray(rotation_cluster_prob_table, [further_left_variable_name, further_right_variable_name, output_variable_name]))
        add_edge_between(variables[further_left_variable_name], factors[factor_name])
        add_edge_between(variables[further_right_variable_name], factors[factor_name])
        add_edge_between(variables[output_variable_name], factors[factor_name])
    end

    number_of_operations["rot"] += 1
    location_execution_counts[input] += 1
end

function xor_with_cluster_shift_factor_graph!(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    input_a::Int64,
    input_b::Int64,
    output::Int64,
    number_of_clusters_shifted::Int64,
    number_of_bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    number_of_operations::Dict{String,Int64},
    precalculated_prob_tables::Dict{String,Array{Float64}},
    round_variables::Set{String},
    round_factors::Set{String},
    run_number::Int64)

    xor_cluster_prob_table = precalculated_prob_tables["xor_cluster"]
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    for i in 1:number_of_clusters
        cur_factor_name = string("f_xor_", number_of_operations["xor"], "_", i, "_", run_number)
        if location_execution_counts[input_a] == 0 && input_a != 13
            input_a_name = string(input_a, "_", location_execution_counts[input_a], "_", i, "_", 1)
        else
            input_a_name = string(input_a, "_", location_execution_counts[input_a], "_", i, "_", run_number)
        end
        if location_execution_counts[input_b] == 0 && input_b != 13
            input_b_name = string(input_b, "_", location_execution_counts[input_b], "_", i, "_", 1)
        else
            input_b_name = string(input_b, "_", location_execution_counts[input_b], "_", i, "_", run_number)
        end
        output_name = string(output, "_", location_execution_counts[output] + 1, "_", ((i - 1 + number_of_clusters + number_of_clusters_shifted) % number_of_clusters) + 1, "_", run_number)
        push!(round_variables, input_a_name)
        push!(round_variables, input_b_name)
        push!(round_variables, output_name)
        push!(round_factors, cur_factor_name)

        factors[cur_factor_name] = XorFactor{AbsVariable}(cur_factor_name)
        variables[output_name] = Variable{AbsFactor}(output_name, number_of_bits_per_cluster)

        # Add the connections between the variables and the factor
        add_edge_between(variables[input_a_name], factors[cur_factor_name])
        add_edge_between(variables[input_b_name], factors[cur_factor_name])
        add_edge_between(variables[output_name], factors[cur_factor_name])
    end
    number_of_operations["xor"] += 1
    location_execution_counts[output] += 1
end

function add_factor_graph!(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    input_a::Int64,
    input_a_version::Int64,
    input_b::Int64,
    input_b_version::Int64,
    output::Int64,
    number_of_bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    number_of_operations::Dict{String,Int64},
    precalculated_prob_tables::Dict{String,Array{Float64}},
    round_variables::Set{String},
    round_factors::Set{String},
    run_number::Int64)

    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    full_add_dist = precalculated_prob_tables["full_add_cluster"]
    add_full_to_output = precalculated_prob_tables["add_full_to_output_cluster"]
    add_full_to_carry = precalculated_prob_tables["add_full_to_carry_cluster"]

    # Could technically make this shared across all the adds but initially for simplicity
    # each will have a seperate one
    initial_carry_var_name = string("add_", number_of_operations["add"], "_carry_0", "_", run_number)
    variables[initial_carry_var_name] = Variable{AbsFactor}(initial_carry_var_name, 1)
    factors["f_"*initial_carry_var_name*"_dist"] = Factor{AbsVariable}("f_" * initial_carry_var_name * "_dist",
        LabelledArray([
                1.0
                0.0
            ], [initial_carry_var_name]))

    for i in 1:number_of_clusters
        full_add_factor_name = string("f_add_", number_of_operations["add"], "_", i, "_", run_number)
        add_carry_out_factor_name = string("f_add_carry_", number_of_operations["add"], "_", i, "_", run_number)
        add_output_factor_name = string("f_add_output_", number_of_operations["add"], "_", i, "_", run_number)

        carry_in_variable_name = string("add_", number_of_operations["add"], "_carry_", i - 1, "_", run_number)
        carry_out_variable_name = string("add_", number_of_operations["add"], "_carry_", i, "_", run_number)
        if input_a_version == 0 && input_a != 13
            input_a_name = string(input_a, "_", input_a_version, "_", i, "_", 1)
        else
            input_a_name = string(input_a, "_", input_a_version, "_", i, "_", run_number)
        end
        if input_b_version == 0 && input_b != 13
            input_b_name = string(input_b, "_", input_b_version, "_", i, "_", 1)
        else
            input_b_name = string(input_b, "_", input_b_version, "_", i, "_", run_number)
        end
        # input_a_name = string(input_a, "_", input_a_version, "_", i, "_", run_number)
        # input_b_name = string(input_b, "_", input_b_version, "_", i, "_", run_number)
        output_name = string(output, "_", location_execution_counts[output] + 1, "_", i, "_", run_number)
        full_add_output_name = string("add_full_out_", number_of_operations["add"], "_", i, "_", run_number)

        push!(round_variables, carry_in_variable_name)
        push!(round_variables, carry_out_variable_name)
        push!(round_variables, input_a_name)
        push!(round_variables, input_b_name)
        push!(round_variables, output_name)
        push!(round_variables, full_add_output_name)
        push!(round_factors, full_add_factor_name)
        push!(round_factors, add_carry_out_factor_name)
        push!(round_factors, add_output_factor_name)

        variables[carry_out_variable_name] = Variable{AbsFactor}(carry_out_variable_name, 1)
        variables[output_name] = Variable{AbsFactor}(output_name, number_of_bits_per_cluster)
        variables[full_add_output_name] = Variable{AbsFactor}(full_add_output_name, number_of_bits_per_cluster + 1)

        factors[full_add_factor_name] = Factor{AbsVariable}(full_add_factor_name, LabelledArray(full_add_dist, [carry_in_variable_name, input_a_name, input_b_name, full_add_output_name]))
        factors[add_carry_out_factor_name] = Factor{AbsVariable}(add_carry_out_factor_name, LabelledArray(add_full_to_carry, [full_add_output_name, carry_out_variable_name]))
        factors[add_output_factor_name] = Factor{AbsVariable}(add_output_factor_name, LabelledArray(add_full_to_output, [full_add_output_name, output_name]))

        add_edge_between(variables[carry_in_variable_name], factors[full_add_factor_name])
        add_edge_between(variables[input_a_name], factors[full_add_factor_name])
        add_edge_between(variables[input_b_name], factors[full_add_factor_name])
        add_edge_between(variables[full_add_output_name], factors[full_add_factor_name])

        add_edge_between(variables[full_add_output_name], factors[add_carry_out_factor_name])
        add_edge_between(variables[carry_out_variable_name], factors[add_carry_out_factor_name])

        add_edge_between(variables[full_add_output_name], factors[add_output_factor_name])
        add_edge_between(variables[output_name], factors[add_output_factor_name])
    end

    #Ensuring that the distribution for the initial carry is the last neighbour for that variable
    add_edge_between(variables[initial_carry_var_name], factors["f_"*initial_carry_var_name*"_dist"])
    number_of_operations["add"] += 1
    location_execution_counts[output] += 1
end

function chacha_quarter_round_factor_graph!(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    number_of_bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    number_of_operations::Dict{String,Int64},
    precalculated_prob_tables::Dict{String,Array{Float64}},
    round_variables::Set{String},
    round_factors::Set{String},
    run_number::Int64)

    first_number_cluster_shifts = Int64(floor(16 / number_of_bits_per_cluster))
    first_number_of_bits_left_to_shift = 16 % number_of_bits_per_cluster

    # println(first_number_cluster_shifts)
    # println(first_number_of_bits_left_to_shift)

    add_factor_graph!(variables, factors, a, location_execution_counts[a], b, location_execution_counts[b], a, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    xor_with_cluster_shift_factor_graph!(variables, factors, d, a, d, first_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    bit_rotation_factor_graph!(variables, factors, d, first_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    # block[a] += block[b]
    # block[d] ⊻= block[a]
    # block[d] = ROTL(block[d], 16)


    second_number_cluster_shifts = Int64(floor(12 / number_of_bits_per_cluster))
    second_number_of_bits_left_to_shift = 12 % number_of_bits_per_cluster

    add_factor_graph!(variables, factors, c, location_execution_counts[c], d, location_execution_counts[d], c, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    xor_with_cluster_shift_factor_graph!(variables, factors, b, c, b, second_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    bit_rotation_factor_graph!(variables, factors, b, second_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    # block[c] += block[d]
    # block[b] ⊻= block[c]
    # block[b] = ROTL(block[b], 12)

    third_number_cluster_shifts = Int64(floor(8 / number_of_bits_per_cluster))
    third_number_of_bits_left_to_shift = 8 % number_of_bits_per_cluster

    add_factor_graph!(variables, factors, a, location_execution_counts[a], b, location_execution_counts[b], a, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    xor_with_cluster_shift_factor_graph!(variables, factors, d, a, d, third_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    bit_rotation_factor_graph!(variables, factors, d, third_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    # block[a] += block[b]
    # block[d] ⊻= block[a]
    # block[d] = ROTL(block[d], 8)

    fourth_number_cluster_shifts = Int64(floor(7 / number_of_bits_per_cluster))
    fourth_number_of_bits_left_to_shift = 7 % number_of_bits_per_cluster

    add_factor_graph!(variables, factors, c, location_execution_counts[c], d, location_execution_counts[d], c, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    xor_with_cluster_shift_factor_graph!(variables, factors, b, c, b, fourth_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    bit_rotation_factor_graph!(variables, factors, b, fourth_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, round_variables, round_factors, run_number)
    # block[c] += block[d]
    # block[b] ⊻= block[c]
    # block[b] = ROTL(block[b], 7)
end

function chacha_factor_graph!(variables::Dict{String, AbsVariable},
    factors::Dict{String,AbsFactor},
    number_of_bits_per_cluster::Int64,
    variables_by_round::Vector{Set{String}},
    factors_by_round::Vector{Set{String}},
    adds_by_round::Vector{Set{Int64}},
    location_execution_counts::Vector{Int64},
    run_number::Int64
)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))

    # Initially add the opening variables
    if !("1_0_1_1" in keys(variables))
        for i in 1:16
            for j in 1:number_of_clusters
                variables[string(i, "_0_", j, "_1")] = Variable{AbsFactor}(string(i, "_0_", j, "_1"), number_of_bits_per_cluster)
            end
        end
    else
        for j in 1:number_of_clusters
            variables[string("13_0_", j, "_", run_number)] = Variable{AbsFactor}(string("13_0_", j, "_", run_number), number_of_bits_per_cluster)
        end
    end
    number_of_operations = Dict("xor" => 0, "add" => 0, "rot" => 0)
    precalculated_prob_tables = Dict("xor_cluster" => make_xor_prob_table(number_of_bits_per_cluster),
        "full_add_cluster" => make_add_including_carry_prob_array(number_of_bits_per_cluster),
        "add_full_to_output_cluster" => take_bottom_bits_prob_array(number_of_bits_per_cluster + 1),
        "add_full_to_carry_cluster" => take_top_bit_prob_array(number_of_bits_per_cluster + 1)
    )

    for j in [16, 12, 8, 7]
        if j % number_of_bits_per_cluster != 0
            precalculated_prob_tables[string("rotation_cluster_", j % number_of_bits_per_cluster)] = make_rotation_prob_table(number_of_bits_per_cluster, j % number_of_bits_per_cluster)
        end
    end

    for i in 1:10
        cur_round_variables = Set{String}()
        cur_round_factors = Set{String}()
        start_add = number_of_operations["add"]
        chacha_quarter_round_factor_graph!(variables, factors, 1, 5, 9, 13, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        chacha_quarter_round_factor_graph!(variables, factors, 2, 6, 10, 14, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        chacha_quarter_round_factor_graph!(variables, factors, 3, 7, 11, 15, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        chacha_quarter_round_factor_graph!(variables, factors, 4, 8, 12, 16, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        variables_by_round[2*i-1] = union(variables_by_round[2*i-1], cur_round_variables)
        factors_by_round[2*i-1] = union(factors_by_round[2*i-1], cur_round_factors)
        adds_by_round[2*i-1] = union(adds_by_round[2*i-1], Set(start_add:(number_of_operations["add"]-1)))
        cur_round_variables = Set{String}()
        cur_round_factors = Set{String}()
        start_add = number_of_operations["add"]
        chacha_quarter_round_factor_graph!(variables, factors, 1, 6, 11, 16, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        chacha_quarter_round_factor_graph!(variables, factors, 2, 7, 12, 13, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        chacha_quarter_round_factor_graph!(variables, factors, 3, 8, 9, 14, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        chacha_quarter_round_factor_graph!(variables, factors, 4, 5, 10, 15, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
        variables_by_round[2*i] = union(variables_by_round[2*i], cur_round_variables)
        factors_by_round[2*i] = union(factors_by_round[2*i], cur_round_factors)
        adds_by_round[2*i] = union(adds_by_round[2*i], Set(start_add:(number_of_operations["add"]-1)))
    end

    # At this stage need to put in the add between the original values and the current value 
    cur_round_variables = Set{String}()
    cur_round_factors = Set{String}()
    start_add = number_of_operations["add"]
    for i in 1:16
        add_factor_graph!(variables, factors, i, 0, i, location_execution_counts[i], i, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables, cur_round_variables, cur_round_factors, run_number)
    end
    variables_by_round[end] = union(variables_by_round[end], cur_round_variables)
    factors_by_round[end] = union(factors_by_round[end], cur_round_factors)
    adds_by_round[end] = union(adds_by_round[end], Set(start_add:(number_of_operations["add"]-1)))
    # This means that all of the variables and factors for the execution of the algorithm have been gone through
    # need to now add in factors for the initial conditions of the ChaCha algorithm this is the first few words being
    # set in the standard
end

function add_starting_constant_values(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    number_of_bits_per_cluster::Int64,
    run_number::Int64)
    set_variable_to_value(variables, factors, string("1_0"), 0x61707865, number_of_bits_per_cluster, run_number)
    set_variable_to_value(variables, factors, string("2_0"), 0x3320646e, number_of_bits_per_cluster, run_number)
    set_variable_to_value(variables, factors, string("3_0"), 0x79622d32, number_of_bits_per_cluster, run_number)
    set_variable_to_value(variables, factors, string("4_0"), 0x6b206574, number_of_bits_per_cluster, run_number)
end

function add_distribution_of_initial_values(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    number_of_bits_per_cluster::Int64,
    key::Vector{UInt32}, nonce::Vector{UInt32}, counter::UInt32,
    run_number::Int64)
    for i in 1:8
        set_variable_to_value(variables, factors, string(i + 4, "_0"), key[i], number_of_bits_per_cluster, run_number)
    end
    set_variable_to_value(variables, factors, string("13_0"), counter, number_of_bits_per_cluster, run_number)
    for i in 1:3
        set_variable_to_value(variables, factors, string(i + 13, "_0"), nonce[i], number_of_bits_per_cluster, run_number)
    end
end

function add_values_of_initial_nonce_and_counter(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    number_of_bits_per_cluster::Int64,
    nonce::Vector{UInt32}, counter::UInt32,
    run_number::Int64)
    set_variable_to_value(variables, factors, string("13_0"), counter, number_of_bits_per_cluster, run_number)
    for i in 1:3
        set_variable_to_value(variables, factors, string(i + 13, "_0"), nonce[i], number_of_bits_per_cluster, run_number)
    end
end

function add_equality_between_keys_and_nonces(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    number_of_bits_per_cluster::Int64,
    number_encryption_runs::Int64,
    factors_names_added::Set{String})

    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    eqaulity_prob_table = make_equality_prob_array(number_of_bits_per_cluster, number_encryption_runs)
    for position_in_state in 5:16
        if position_in_state != 13
            for i in 1:number_of_clusters
                factor_name = string("f_", position_in_state, "_0_", i, "_equality")
                variable_names = [string(position_in_state, "_0_", i, "_", run_number) for run_number in 1:number_encryption_runs]
                factors[factor_name] = Factor{AbsVariable}(factor_name, LabelledArray(eqaulity_prob_table, variable_names))
                for var_name in variable_names
                    add_edge_between(variables[var_name], factors[factor_name])
                end
                push!(factors_names_added, factor_name)
            end
        end
    end
end

function add_adds_between_counters(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    number_of_bits_per_cluster::Int64,
    number_encryption_runs::Int64,
    factors_names_added::Set{String},
    variable_names_add::Set{String})

    # Make an add between the counters (the 16th item of the state) of the consecutive runs of the algorithm because
    # we know that it is incremented by 1 between each run
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    full_add_dist = make_add_including_carry_prob_array(number_of_bits_per_cluster)
    add_full_to_output = take_bottom_bits_prob_array(number_of_bits_per_cluster + 1)
    add_full_to_carry = take_top_bit_prob_array(number_of_bits_per_cluster + 1)

    zero_dist_table = zeros(1 << number_of_bits_per_cluster)
    zero_dist_table[1] = 1.0

    for run_number in 1:(number_encryption_runs-1)
        # Need to create a variable which has a value equal to 1 so that it can be added to the value
        # but it is slightly simpler to just set the carry bit and add a value of zero and that should 
        # lead to an identical effect

        initial_carry_var_name = string("add_between_runs_", run_number, "_carry_0")
        variables[initial_carry_var_name] = Variable{AbsFactor}(initial_carry_var_name, 1)
        factors["f_"*initial_carry_var_name*"_dist"] = Factor{AbsVariable}("f_" * initial_carry_var_name * "_dist",
            LabelledArray([
                    0.0
                    1.0
                ], [initial_carry_var_name]))
        add_edge_between(variables[initial_carry_var_name], factors["f_"*initial_carry_var_name*"_dist"])

        for i in 1:number_of_clusters
            full_add_factor_name = string("f_add_between_runs_", run_number, "_", i)
            add_carry_out_factor_name = string("f_add_between_runs_carry_", run_number, "_", i)
            add_output_factor_name = string("f_add_between_runs_output_", run_number, "_", i)

            carry_in_variable_name = string("add_between_runs_", run_number, "_carry_", i - 1)
            carry_out_variable_name = string("add_between_runs_", run_number, "_carry_", i)
            input_a_name = string("13_0_", i, "_", run_number)

            input_b_name = string("zero_", i, "_", run_number)
            input_b_dist_name = string("f_", input_b_name, "_dist")
            # Don't need to add these to the set of values to be evaluated on every iteration because it will not change
            # due to just being a constant zero
            variables[input_b_name] = Variable{AbsFactor}(input_b_name, number_of_bits_per_cluster)
            factors[input_b_dist_name] = Factor{AbsVariable}(input_b_dist_name, LabelledArray(zero_dist_table, [input_b_name]))
            add_edge_between(variables[input_b_name], factors[input_b_dist_name])

            output_name = string("13_0_", i, "_", run_number + 1)
            full_add_output_name = string("add_between_runs_", i, "_", run_number)

            push!(variable_names_add, carry_in_variable_name)
            push!(variable_names_add, carry_out_variable_name)
            push!(variable_names_add, full_add_output_name)
            push!(factors_names_added, full_add_factor_name)
            push!(factors_names_added, add_carry_out_factor_name)
            push!(factors_names_added, add_output_factor_name)

            variables[carry_out_variable_name] = Variable{AbsFactor}(carry_out_variable_name, 1)
            variables[full_add_output_name] = Variable{AbsFactor}(full_add_output_name, number_of_bits_per_cluster + 1)

            factors[full_add_factor_name] = Factor{AbsVariable}(full_add_factor_name, LabelledArray(full_add_dist, [carry_in_variable_name, input_a_name, input_b_name, full_add_output_name]))
            factors[add_carry_out_factor_name] = Factor{AbsVariable}(add_carry_out_factor_name, LabelledArray(add_full_to_carry, [full_add_output_name, carry_out_variable_name]))
            factors[add_output_factor_name] = Factor{AbsVariable}(add_output_factor_name, LabelledArray(add_full_to_output, [full_add_output_name, output_name]))

            add_edge_between(variables[carry_in_variable_name], factors[full_add_factor_name])
            add_edge_between(variables[input_a_name], factors[full_add_factor_name])
            add_edge_between(variables[input_b_name], factors[full_add_factor_name])
            add_edge_between(variables[full_add_output_name], factors[full_add_factor_name])

            add_edge_between(variables[full_add_output_name], factors[add_carry_out_factor_name])
            add_edge_between(variables[carry_out_variable_name], factors[add_carry_out_factor_name])

            add_edge_between(variables[full_add_output_name], factors[add_output_factor_name])
            add_edge_between(variables[output_name], factors[add_output_factor_name])
        end
    end
end

function belief_propagate_through_add(variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    bits_per_cluster::Int64,
    add_number::Int64,
    run_number::Int64)

    number_of_clusters = Int64(ceil(32 / bits_per_cluster))
    does_not_contain_add(x::String) = !occursin("add", x)
    inputs = filter(does_not_contain_add, [variable.name for variable in factors[string("f_add_", add_number, "_1", "_", run_number)].neighbours])
    output = filter(does_not_contain_add, [variable.name for variable in factors[string("f_add_output_", add_number, "_1", "_", run_number)].neighbours])
    # Initially pass messages from all of the inputs
    for i in 1:number_of_clusters
        for cur_input in inputs
            variable_to_factor_messages(variables[string(cur_input[begin:end-3], i, "_", run_number)])
        end
    end
    # Then pass messages from all the outputs except the last
    for i in 1:(number_of_clusters-1)
        for cur_output in output
            variable_to_factor_messages(variables[string(cur_output[begin:end-3], i, "_", run_number)])
        end
    end
    # Also pass a message from the final carry
    variable_to_factor_messages(variables[string("add_", add_number, "_carry_", number_of_clusters, "_", run_number)])

    # Now do the heart of the process of going through the tree passing messages
    for i in 1:(number_of_clusters-1)
        variable_to_factor_messages(variables[string("add_", add_number, "_carry_", i - 1, "_", run_number)])
        factor_to_variable_messages(factors[string("f_add_output_", add_number, "_", i, "_", run_number)])
        factor_to_variable_messages(factors[string("f_add_", add_number, "_", i, "_", run_number)])
        variable_to_factor_messages(variables[string("add_full_out_", add_number, "_", i, "_", run_number)])
        factor_to_variable_messages(factors[string("f_add_carry_", add_number, "_", i, "_", run_number)])
    end

    # Perform the passing messages for the final cluster
    variable_to_factor_messages(variables[string("add_", add_number, "_carry_", number_of_clusters - 1, "_", run_number)])
    factor_to_variable_messages(factors[string("f_add_", add_number, "_", number_of_clusters, "_", run_number)])
    factor_to_variable_messages(factors[string("f_add_carry_", add_number, "_", number_of_clusters, "_", run_number)])
    variable_to_factor_messages(variables[string("add_full_out_", add_number, "_", number_of_clusters, "_", run_number)])
    factor_to_variable_messages(factors[string("f_add_output_", add_number, "_", number_of_clusters, "_", run_number)])

    for cur_output in output
        variable_to_factor_messages(variables[string(cur_output[begin:end-3], number_of_clusters, "_", run_number)])
    end

    # Now perform the same message passing in reverse which should give the correct marginalisation through the add
    factor_to_variable_messages(factors[string("f_add_output_", add_number, "_", number_of_clusters, "_", run_number)])
    variable_to_factor_messages(variables[string("add_full_out_", add_number, "_", number_of_clusters, "_", run_number)])
    factor_to_variable_messages(factors[string("f_add_carry_", add_number, "_", number_of_clusters, "_", run_number)])
    factor_to_variable_messages(factors[string("f_add_", add_number, "_", number_of_clusters, "_", run_number)])
    variable_to_factor_messages(variables[string("add_", add_number, "_carry_", number_of_clusters - 1, "_", run_number)])

    # Now do the heart of the process of going through the tree passing messages
    for i in (number_of_clusters-1):1
        factor_to_variable_messages(factors[string("f_add_carry_", add_number, "_", i, "_", run_number)])
        variable_to_factor_messages(variables[string("add_full_out_", add_number, "_", i, "_", run_number)])
        factor_to_variable_messages(factors[string("f_add_", add_number, "_", i, "_", run_number)])
        factor_to_variable_messages(factors[string("f_add_output_", add_number, "_", i, "_", run_number)])
        variable_to_factor_messages(variables[string("add_", add_number, "_carry_", i - 1, "_", run_number)])
    end

    variable_to_factor_messages(variables[string("add_", add_number, "_carry_", number_of_clusters, "_", run_number)])
    for i in (number_of_clusters-1):1
        for cur_output in output
            variable_to_factor_messages(variables[string(cur_output[begin:end-3], i, "_", run_number)])
        end
    end
    for i in number_of_clusters:1
        for cur_input in inputs
            variable_to_factor_messages(variables[string(cur_input[begin:end-3], i, "_", run_number)])
        end
    end
end