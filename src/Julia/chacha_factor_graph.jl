include("messages.jl")

# This creates an entire factor graph for the ChaCha algorithm

# It works by having a numbering for the different locations in ChaCha algorithm as 
# 1  2  3  4 
# 5  6  7  8 
# 9  10 11 12
# 13 14 15 16

# Then going forward the variables will be labeled as {location}_{number_of_times_used}_{cluster_number}, 
# with the additional variables introduced by variables but not native in the algorithm (e.g. carry bits, temporary parts in adds)
# these will also be labeled by refering to the number of times they have been used alongside things such as location

function make_rotation_prob_table(number_of_bits_in_cluster::Int64, number_of_bits_to_shift::Int64)
    output = zeros(1 << number_of_bits_in_cluster, 1 << number_of_bits_in_cluster, 1 << number_of_bits_in_cluster)
    number_of_bits_from_right = number_of_bits_to_shift
    number_of_bits_from_left = number_of_bits_in_cluster - number_of_bits_from_right
    for input_left in 1:(1 << number_of_bits_in_cluster)
        for input_right in 1:(1 << number_of_bits_in_cluster)
            output[input_left, input_right, ((((input_left - 1) & ((1 << number_of_bits_from_left) - 1)) << number_of_bits_from_right)|((((input_right >> number_of_bits_from_left) - 1) & ((1 << number_of_bits_from_right) - 1)))) + 1] = 1.
        end
    end
    return output
end

function make_xor_prob_table(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 1 << num_of_bits, 1 << num_of_bits)
    for input_a in 1:(1 << num_of_bits)
        for input_b in 1:(1 << num_of_bits)
            output[input_a, input_b, ((input_a - 1) ⊻ (input_b - 1)) + 1] = 1.
        end
    end
    return output
end


function make_add_including_carry_prob_array(num_of_bits::Int64)
    output = zeros(2, 1 << num_of_bits, 1 << num_of_bits, 1 << (num_of_bits + 1))
    for carry_in in 1:2
        for input_a in 1:(1 << num_of_bits)
            for input_b in 1:(1 << num_of_bits)
                output[carry_in, input_a, input_b, input_a + input_b + carry_in - 2] = 1.
            end
        end
    end
    return output
end

function take_top_bit_prob_array(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 2)
    for val_in in 1:(1 << num_of_bits)
        for output_val in 1:2
            if (output_val == 2 && val_in > 1 << (num_of_bits - 1)) || (output_val == 1 && val_in <= 1 << (num_of_bits - 1))
                output[val_in, output_val] = 1.
            end
        end
    end
    return output
end

function take_bottom_bits_prob_array(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 1 << (num_of_bits - 1))
    for val_in in 1:1 << num_of_bits
        output[val_in, ((val_in - 1) % ((1 << (num_of_bits - 1)))) + 1] = 1.
    end
    return output
end

function bit_rotation_factor_graph!(variables,
    factors,
    input,
    bits_to_rotate_by,
    number_of_bits_per_cluster,
    location_execution_counts,
    number_of_operations,
    precalculated_prob_tables
)
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
        further_left_variable_name = string(input, "_", location_execution_counts[input], "_", i)
        further_right_variable_name = string(input, "_", location_execution_counts[input], "_", i - 1 <= 0 ? i - 1 + number_of_clusters : i - 1)
        output_variable_name = string(input, "_", location_execution_counts[input] + 1, "_", i)
        factor_name = string("f_rot_", number_of_operations["rot"], "_", i)
        
        variables[output_variable_name] = Variable(output_variable_name)
        factors[factor_name] = Factor(factor_name, LabelledArray(rotation_cluster_prob_table, [further_left_variable_name, further_right_variable_name, output_variable_name]))
        add_edge_between(variables[further_left_variable_name], factors[factor_name])
        add_edge_between(variables[further_right_variable_name], factors[factor_name])
        add_edge_between(variables[output_variable_name], factors[factor_name])
    end
    
    number_of_operations["rot"] += 1
    location_execution_counts[input] += 1
end

function xor_with_cluster_shift_factor_graph!(variables,
    factors,
    input_a,
    input_b,
    output,
    number_of_clusters_shifted,
    number_of_bits_per_cluster,
    location_execution_counts,
    number_of_operations,
    precalculated_prob_tables)

    xor_cluster_prob_table = precalculated_prob_tables["xor_cluster"]
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    for i in 1:number_of_clusters
        cur_factor_name = string("f_xor_", number_of_operations["xor"], "_", i)
        input_a_name = string(input_a, "_", location_execution_counts[input_a], "_", i)
        input_b_name = string(input_b, "_", location_execution_counts[input_b], "_", i)
        output_name = string(output, "_", location_execution_counts[output] + 1, "_", ((i - 1 + number_of_clusters + number_of_clusters_shifted) % number_of_clusters) + 1)
        factors[cur_factor_name] = Factor(cur_factor_name, LabelledArray(xor_cluster_prob_table,
        [input_a_name, input_b_name, output_name]))
        variables[output_name] = Variable(output_name)

        # Add the connections between the variables and the factor
        add_edge_between(variables[input_a_name], factors[cur_factor_name])
        add_edge_between(variables[input_b_name], factors[cur_factor_name])
        add_edge_between(variables[output_name], factors[cur_factor_name])
    end
    number_of_operations["xor"] += 1
    location_execution_counts[output] += 1
end

function add_factor_graph!(variables,
    factors,
    input_a,
    input_a_version,
    input_b,
    input_b_version,
    output,
    number_of_bits_per_cluster,
    location_execution_counts,
    number_of_operations,
    precalculated_prob_tables)

    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    full_add_dist = precalculated_prob_tables["full_add_cluster"]
    add_full_to_output = precalculated_prob_tables["add_full_to_output_cluster"]
    add_full_to_carry = precalculated_prob_tables["add_full_to_carry_cluster"]

    # Could technically make this shared across all the adds but initially for simplicity
    # each will have a seperate one
    initial_carry_var_name = string("add_", number_of_operations["add"], "_carry_0")
    variables[initial_carry_var_name] = Variable(initial_carry_var_name)
    factors["f_" * initial_carry_var_name * "_dist"] = Factor("f_" * initial_carry_var_name * "_dist",
        LabelledArray([
            1.
            0.
        ], [initial_carry_var_name]))
    add_edge_between(variables[initial_carry_var_name], factors["f_" * initial_carry_var_name * "_dist"])

    for i in 1:number_of_clusters
        full_add_factor_name = string("f_add_", number_of_operations["add"], "_", i)
        add_carry_out_factor_name = string("f_add_carry_", number_of_operations["add"], "_", i)
        add_output_factor_name = string("f_add_output_", number_of_operations["add"], "_", i)

        carry_in_variable_name = string("add_", number_of_operations["add"], "_carry_", i - 1)
        carry_out_variable_name = string("add_", number_of_operations["add"], "_carry_", i)
        input_a_name = string(input_a, "_", input_a_version, "_", i)
        input_b_name = string(input_b, "_", input_b_version, "_", i)
        output_name = string(output, "_", location_execution_counts[output] + 1, "_", i)
        full_add_output_name = string("add_full_out_", number_of_operations["add"], "_", i)

        variables[carry_out_variable_name] = Variable(carry_out_variable_name)
        variables[output_name] = Variable(output_name)
        variables[full_add_output_name] = Variable(full_add_output_name)

        factors[full_add_factor_name] = Factor(full_add_factor_name, LabelledArray(full_add_dist, [carry_in_variable_name, input_a_name, input_b_name, output_name]))
        factors[add_carry_out_factor_name] = Factor(add_carry_out_factor_name, LabelledArray(add_full_to_carry, [full_add_output_name, carry_out_variable_name]))
        factors[add_output_factor_name] = Factor(add_output_factor_name, LabelledArray(add_full_to_output, [full_add_output_name, output_name]))

        add_edge_between(variables[carry_in_variable_name], factors[full_add_factor_name])
        add_edge_between(variables[input_a_name], factors[full_add_factor_name])
        add_edge_between(variables[input_b_name], factors[full_add_factor_name])
        add_edge_between(variables[full_add_output_name], factors[full_add_factor_name])

        add_edge_between(variables[full_add_output_name], factors[add_carry_out_factor_name])
        add_edge_between(variables[carry_out_variable_name], factors[add_carry_out_factor_name])

        add_edge_between(variables[full_add_output_name], factors[add_output_factor_name])
        add_edge_between(variables[output_name], factors[add_output_factor_name])
    end

    number_of_operations["add"] += 1
    location_execution_counts[output] += 1
end

function chacha_quarter_round_factor_graph!(variables,
    factors,
    a,
    b,
    c,
    d,
    number_of_bits_per_cluster,
    location_execution_counts,
    number_of_operations,
    precalculated_prob_tables)
    
    first_number_cluster_shifts = Int64(floor(16 / number_of_bits_per_cluster))
    first_number_of_bits_left_to_shift = 16 % number_of_bits_per_cluster

    # println(first_number_cluster_shifts)
    # println(first_number_of_bits_left_to_shift)

    add_factor_graph!(variables, factors, a, location_execution_counts[a], b, location_execution_counts[b], a, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, d, a, d, first_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, d, first_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[a] += block[b]
    # block[d] ⊻= block[a]
    # block[d] = ROTL(block[d], 16)


    second_number_cluster_shifts = Int64(floor(12 / number_of_bits_per_cluster))
    second_number_of_bits_left_to_shift = 12 % number_of_bits_per_cluster
    
    add_factor_graph!(variables, factors, c, location_execution_counts[c], d, location_execution_counts[d], c, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, b, c, b, second_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, b, second_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[c] += block[d]
    # block[b] ⊻= block[c]
    # block[b] = ROTL(block[b], 12)

    third_number_cluster_shifts = Int64(floor(8 / number_of_bits_per_cluster))
    third_number_of_bits_left_to_shift = 8 % number_of_bits_per_cluster
    
    add_factor_graph!(variables, factors, a, location_execution_counts[a], b, location_execution_counts[b], a, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, d, a, d, third_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, d, third_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[a] += block[b]
    # block[d] ⊻= block[a]
    # block[d] = ROTL(block[d], 8)

    fourth_number_cluster_shifts = Int64(floor(7 / number_of_bits_per_cluster))
    fourth_number_of_bits_left_to_shift = 7 % number_of_bits_per_cluster
    
    add_factor_graph!(variables, factors, c, location_execution_counts[c], d, location_execution_counts[d], c, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, b, c, b, fourth_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, b, fourth_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[c] += block[d]
    # block[b] ⊻= block[c]
    # block[b] = ROTL(block[b], 7)
end

function chacha_factor_graph!(variables,
    factors,
    number_of_bits_per_cluster)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    
    # Initially add the opening variables
    for i in 1:16
        for j in 1:number_of_clusters
            variables[string(i, "_0_", j)] = Variable(string(i, "_0_", j))
        end
    end

    location_execution_counts = zeros(Int64, 16)
    number_of_operations = Dict("xor" => 0, "add" => 0, "rot" => 0)
    precalculated_prob_tables = Dict("xor_cluster" => make_xor_prob_table(number_of_bits_per_cluster),
    "full_add_cluster" => make_add_including_carry_prob_array(number_of_bits_per_cluster),
    "add_full_to_output_cluster" =>  take_bottom_bits_prob_array(number_of_bits_per_cluster + 1),
    "add_full_to_carry_cluster" => take_top_bit_prob_array(number_of_bits_per_cluster + 1)
    )

    for j in [16, 12, 8, 7]
        if j % number_of_bits_per_cluster != 0
            precalculated_prob_tables[string("rotation_cluster_", j % number_of_bits_per_cluster)] = make_rotation_prob_table(number_of_bits_per_cluster, j % number_of_bits_per_cluster)
        end
    end

    for i in 1:10
        chacha_quarter_round_factor_graph!(variables, factors, 1, 5, 9, 13, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        chacha_quarter_round_factor_graph!(variables, factors, 2, 6, 10, 14, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        chacha_quarter_round_factor_graph!(variables, factors, 3, 7, 11, 15, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        chacha_quarter_round_factor_graph!(variables, factors, 4, 8, 12, 16, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        
        chacha_quarter_round_factor_graph!(variables, factors, 1, 6, 11, 16, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        chacha_quarter_round_factor_graph!(variables, factors, 2, 7, 12, 13, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        chacha_quarter_round_factor_graph!(variables, factors, 3, 8, 9, 14, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
        chacha_quarter_round_factor_graph!(variables, factors, 4, 5, 10, 15, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    end

    # At this stage need to put in the add between the original values and the current value 
    for i in 1:16
        add_factor_graph!(variables, factors, i, 0, i, location_execution_counts[i], i, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    end

    # This means that all of the variables and factors for the execution of the algorithm have been gone through
    # need to now add in factors for the initial conditions of the ChaCha algorithm this is the first few words being
    # set in the standard
end

number_of_bits = 2
variables = Dict()
factors = Dict()
chacha_factor_graph!(variables, factors, number_of_bits)

for (i,j) in factors
    println(i)
    factor_to_variable_messages(j)
end

for (i,j) in variables
    variable_to_factor_messages(j)
end
