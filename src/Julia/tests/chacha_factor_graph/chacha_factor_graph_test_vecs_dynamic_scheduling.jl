using DataStructures
using Test

include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../belief_propagation/dynamic_message_scheduling.jl")

function test_quarter_round()
    bits_per_cluster = 2
    number_of_clusters = Int64(ceil(32 / bits_per_cluster))
    initial_state = [0x879531e0, 0xc5ecf37d, 0x516461b1, 0xc9a62f8a, 0x44c20ef3, 0x3390af7f, 0xd9fc690b, 0x2a5f714c, 0x53372767, 0xb00a5631, 0x974c541a, 0x359e9963, 0x5c971061, 0x3d631689, 0x2098d9d6, 0x91dbd320]
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    for i in 1:16
        for j in 1:number_of_clusters
            variables[string(i, "_0_", j, "_1")] = Variable{AbsFactor}(string(i, "_0_", j, "_1"), bits_per_cluster)
        end
    end
    location_execution_counts = zeros(Int64, 16)
    number_of_operations = Dict("xor" => 0, "add" => 0, "rot" => 0)
    precalculated_prob_tables = Dict("xor_cluster" => make_xor_prob_table(bits_per_cluster),
    "full_add_cluster" => make_add_including_carry_prob_array(bits_per_cluster),
    "add_full_to_output_cluster" =>  take_bottom_bits_prob_array(bits_per_cluster + 1),
    "add_full_to_carry_cluster" => take_top_bit_prob_array(bits_per_cluster + 1)
    )

    for j in [16, 12, 8, 7]
        if j % bits_per_cluster != 0
            precalculated_prob_tables[string("rotation_cluster_", j % bits_per_cluster)] = make_rotation_prob_table(bits_per_cluster, j % bits_per_cluster)
        end
    end

    round_variables = Set{String}()
    round_factors = Set{String}()
    chacha_quarter_round_factor_graph!(variables, factors, 3, 8, 9, 14, bits_per_cluster, location_execution_counts, number_of_operations, precalculated_prob_tables,
    round_variables, round_factors, 1)
    for i in 1:16
        set_variable_to_value(variables, factors, string(i, "_0"), initial_state[i], bits_per_cluster, 1)
    end

    # Iterate the message passing lots of times to pass data through the graph to hopefully have coverged to the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 5000)

    expected_state = [0x879531e0, 0xc5ecf37d, 0xbdb886dc, 0xc9a62f8a, 0x44c20ef3, 0x3390af7f, 0xd9fc690b, 0xcfacafd2, 0xe46bea80, 0xb00a5631, 0x974c541a, 0x359e9963, 0x5c971061, 0xccc07c79, 0x2098d9d6, 0x91dbd320]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end

function test_vec_sec_2_full_encryption()
    key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
    nonce = [0x09000000, 0x4a000000, 0x00000000]
    counter::UInt32 = 1

    bits_per_cluster = 2
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)

    chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)

    add_starting_constant_values(variables, factors, bits_per_cluster, 1)
    add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, 1)

    # Need to pass messages through the factor graph to reach the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)

    expected_state = [0xe4e7f110, 0x15593bd1, 0x1fdd0f50, 0xc47120a3, 0xc7f4d1c7, 0x0368c033, 0x9aaa2204, 0x4e6cd4c3, 0x466482d2, 0x09aa9f07, 0x05d7c214, 0xa2028bd9, 0xd19c12b5, 0xb94e16de, 0xe883d0cb, 0x4e3c50a2]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end

function test_vec_1_full_encryption()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 0

    bits_per_cluster = 1
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)

    chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)

    add_starting_constant_values(variables, factors, bits_per_cluster, 1)
    add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, 1)

    # Need to pass messages through the factor graph to reach the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)

    expected_state = [0xade0b876, 0x903df1a0, 0xe56a5d40, 0x28bd8653, 0xb819d2bd, 0x1aed8da0, 0xccef36a8, 0xc70d778b, 0x7c5941da, 0x8d485751, 0x3fe02477, 0x374ad8b8, 0xf4b8436a, 0x1ca11815, 0x69b687c3, 0x8665eeb2]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end


function test_vec_2_full_encryption()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 1

    bits_per_cluster = 4
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)

    chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)

    add_starting_constant_values(variables, factors, bits_per_cluster, 1)
    add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, 1)

    # Need to pass messages through the factor graph to reach the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)

    expected_state = [0xbee7079f, 0x7a385155, 0x7c97ba98, 0x0d082d73, 0xa0290fcb, 0x6965e348, 0x3e53c612, 0xed7aee32, 0x7621b729, 0x434ee69c, 0xb03371d5, 0xd539d874, 0x281fed31, 0x45fb0a51, 0x1f0ae1ac, 0x6f4d794b]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end

function test_vec_3_full_encryption()
    key = zeros(UInt32, 8)
    key[8] = 0x01000000
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 1

    bits_per_cluster = 2
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)

    chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)

    add_starting_constant_values(variables, factors, bits_per_cluster, 1)
    add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, 1)

    # Need to pass messages through the factor graph to reach the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)

    expected_state = [0x2452eb3a, 0x9249f8ec, 0x8d829d9b, 0xddd4ceb1, 0xe8252083, 0x60818b01, 0xf38422b8, 0x5aaa49c9, 0xbb00ca8e, 0xda3ba7b4, 0xc4b592d1, 0xfdf2732f, 0x4436274e, 0x2561b3c8, 0xebdd4aa6, 0xa0136c00]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end

function test_vec_4_full_encryption()
    key = zeros(UInt32, 8)
    key[1] = 0x0000ff00
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 2

    bits_per_cluster = 2
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)

    chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)

    add_starting_constant_values(variables, factors, bits_per_cluster, 1)
    add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, 1)

    # Need to pass messages through the factor graph to reach the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)

    expected_state = [0xfb4dd572, 0x4bc42ef1, 0xdf922636, 0x327f1394, 0xa78dea8f, 0x5e269039, 0xa1bebbc1, 0xcaf09aae, 0xa25ab213, 0x48a6b46c, 0x1b9d9bcb, 0x092c5be6, 0x546ca624, 0x1bec45d5, 0x87f47473, 0x96f0992e]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end

function test_vec_5_full_encryption()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    nonce[3] = 0x02000000
    counter::UInt32 = 0

    bits_per_cluster = 2
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)

    chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)

    add_starting_constant_values(variables, factors, bits_per_cluster, 1)
    add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, 1)

    # Need to pass messages through the factor graph to reach the correct values
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)

    expected_state = [0x374dc6c2, 0x3736d58c, 0xb904e24a, 0xcd3f93ef, 0x88228b1a, 0x96a4dfb3, 0x5b76ab72, 0xc727ee54, 0x0e0e978a, 0xf3145c95, 0x1b748ea8, 0xf786c297, 0x99c28f5f, 0x628314e8, 0x398a19fa, 0x6ded1b53]
    actual_most_likely_state = [read_most_likely_value_from_variable(variables, string(i,"_", location_execution_counts[i]), bits_per_cluster, 1) for i in 1:16]
    println(actual_most_likely_state)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_state == actual_most_likely_state
end

# I think that a good strategy for ordering the passing of information around would be 
# to have a priority queue, where the variables with the highest change in entropy are first,
# and then select the highest and do the message passing for it and its neighbours recomputing
# the changes in entropy accordingly updating the priority queue allowing for a dynamic moving of
# information around the graph (may potentially want some additional structure for triggering a better
# passing through the and trees but may not be required.  Then can take the total entropy of the graph
# every few hundred/thousand updates to see how it is fairing (which could now just be a sum rather than)
# needing to calculate it on every iteration

test_quarter_round()
test_vec_sec_2_full_encryption()
test_vec_1_full_encryption()
test_vec_2_full_encryption()
test_vec_3_full_encryption()
test_vec_4_full_encryption()
test_vec_5_full_encryption()