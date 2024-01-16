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

function bit_rotation_factor_graph!(variables,
    factors,
    input,
    bits_to_rotate_by,
    number_of_bits_per_cluster,
    location_execution_counts,
    number_of_operations,
    precalculated_prob_tables
)
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

end

function add_factor_graph!(variables,
    factors,
    input_a,
    input_b,
    output,
    number_of_bits_per_cluster,
    location_execution_counts,
    number_of_operations,
    precalculated_prob_tables)

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
    
    add_factor_graph!(variables, factors, a, b, a, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, d, a, d, first_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, d, first_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[a] += block[b]
    # block[d] ⊻= block[a]
    # block[d] = ROTL(block[d], 16)


    second_number_cluster_shifts = Int64(floor(12 / number_of_bits_per_cluster))
    second_number_of_bits_left_to_shift = 12 % number_of_bits_per_cluster
    
    add_factor_graph!(variables, factors, c, d, c, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, b, c, b, second_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, b, second_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[c] += block[d]
    # block[b] ⊻= block[c]
    # block[b] = ROTL(block[b], 12)

    third_number_cluster_shifts = Int64(floor(8 / number_of_bits_per_cluster))
    third_number_of_bits_left_to_shift = 8 % number_of_bits_per_cluster
    
    add_factor_graph!(variables, factors, a, b, a, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, d, a, d, third_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, d, third_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[a] += block[b]
    # block[d] ⊻= block[a]
    # block[d] = ROTL(block[d], 8)

    fourth_number_cluster_shifts = Int64(floor(7 / number_of_bits_per_cluster))
    fourth_number_of_bits_left_to_shift = 7 % number_of_bits_per_cluster
    
    add_factor_graph!(variables, factors, c, d, c, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    xor_with_cluster_shift_factor_graph!(variables, factors, b, c, b, fourth_number_cluster_shifts, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    bit_rotation_factor_graph!(variables, factors, b, fourth_number_of_bits_left_to_shift, number_of_bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables)
    # block[c] += block[d]
    # block[b] ⊻= block[c]
    # block[b] = ROTL(block[b], 7)
end

function chacha_factor_graph!(variables, factors, number_of_bits_per_cluster)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    
    # Initially add the opening variables
    for i in 1:16
        for j in 1:number_of_clusters
            variables[string(i, "_0_", j)] = Variable(string(i, "_0_", j))
        end
    end

    location_execution_counts = zeros(Int64, 16)
    number_of_operations = Dict("xor" => 0, "add" => 0, "rot" => 0)
    precalculated_prob_tables = Dict()

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


    # This means that all of the variables and factors for the execution of the algorithm have been gone through
    # need to now add in factors for the initial conditions of the ChaCha algorithm this is the first few words being
    # set in the standard
end