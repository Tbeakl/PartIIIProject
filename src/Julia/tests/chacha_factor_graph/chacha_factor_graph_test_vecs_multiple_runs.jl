using Test
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../encryption/chacha.jl")

function belief_propogate_through_graph_forwards(variables::Dict{String,Variable{Factor}},
    factors::Dict{String,Factor{Variable}},
    variables_by_round::Vector{Set{String}},
    factors_by_round::Vector{Set{String}},
    adds_by_round::Vector{Set{Int64}},
    bits_per_cluster::Int64,
    number_of_encryption_runs::Int64
)
    all_variables = [keys(variables)...]
    update_all_entropies(variables, all_variables)
    for (i, j) in factors
        # println(i)
        factor_to_variable_messages(j)
    end
    for (i, j) in variables
        # println(i)
        variable_to_factor_messages(j)
    end
    update_all_entropies(variables, all_variables)
    
    for round_iter in 1:50
        println(round_iter)
        for factor_name in factors_by_round[1]
            # println(factor_name)
            factor_to_variable_messages(factors[factor_name])
        end
        for variable_name in variables_by_round[1]
            # println(variable_name)
            variable_to_factor_messages(variables[variable_name])
        end
    end

    for passes_through_graph in 1:1
        println("On iteration ", passes_through_graph)
        for i in 1:length(variables_by_round)
            println("Round ", i)
            for round_iter in 1:10
                println(round_iter)
                for factor_name in factors_by_round[i]
                    # println(factor_name)
                    factor_to_variable_messages(factors[factor_name])
                end
                for variable_name in variables_by_round[i]
                    # println(variable_name)
                    variable_to_factor_messages(variables[variable_name])
                end
                for add_num in adds_by_round[i]
                    for encryption_run in 1:number_of_encryption_runs
                        belief_propagate_through_add(variables, factors, bits_per_cluster, add_num, encryption_run)
                    end
                end
            end
            update_all_entropies(variables, all_variables)
            println("Total entropy after round ", total_entropy_of_graph(variables))
        end
        println("Total entropy after pass ", total_entropy_of_graph(variables))
    end
end

function test_vec_sec_2_full_encryption()
    key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
    nonce = [0x09000000, 0x4a000000, 0x00000000]
    counter::UInt32 = 1

    bits_per_cluster = 2
    number_encryption_runs = 3
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    expected_states::Vector{Vector{UInt32}} = []

    for i in 1:number_encryption_runs
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, i)
        if i == 1
            add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, i)
        end
        add_starting_constant_values(variables, factors, bits_per_cluster, i)
        push!(expected_states, encrypt(key, nonce, counter))
        counter = counter + 0x1
    end

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()
    # println(number_encryption_runs)
    add_equality_between_keys_and_nonces(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors)
    add_adds_between_counters(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors, additional_variables)

    variables_by_round[1] = union(variables_by_round[1], additional_variables)
    factors_by_round[1] = union(factors_by_round[1], additional_factors)

    # Need to pass messages through the factor graph to reach the correct values
    belief_propogate_through_graph_forwards(variables, factors, variables_by_round, factors_by_round, adds_by_round, bits_per_cluster, number_encryption_runs)

    actual_most_likely_states = [[read_most_likely_value_from_variable(variables, string(i, "_", location_execution_counts[i]), bits_per_cluster, encryption_run) for i in 1:16] for encryption_run in 1:number_encryption_runs]
    println(actual_most_likely_states)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_states == actual_most_likely_states
end

function test_vec_1_full_encryption()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 0

    bits_per_cluster = 1
    number_encryption_runs = 4
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    expected_states::Vector{Vector{UInt32}} = []

    for i in 1:number_encryption_runs
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, i)
        if i == 1
            add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, i)
        end
        add_starting_constant_values(variables, factors, bits_per_cluster, i)
        push!(expected_states, encrypt(key, nonce, counter))
        counter = counter + 0x1
    end

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()
    # println(number_encryption_runs)
    add_equality_between_keys_and_nonces(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors)
    add_adds_between_counters(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors, additional_variables)

    variables_by_round[1] = union(variables_by_round[1], additional_variables)
    factors_by_round[1] = union(factors_by_round[1], additional_factors)

    # Need to pass messages through the factor graph to reach the correct values
    belief_propogate_through_graph_forwards(variables, factors, variables_by_round, factors_by_round, adds_by_round, bits_per_cluster, number_encryption_runs)

    actual_most_likely_states = [[read_most_likely_value_from_variable(variables, string(i, "_", location_execution_counts[i]), bits_per_cluster, encryption_run) for i in 1:16] for encryption_run in 1:number_encryption_runs]
    println(actual_most_likely_states)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_states == actual_most_likely_states
end

function test_vec_2_full_encryption()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 1

    bits_per_cluster = 4
    number_encryption_runs = 2
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    expected_states::Vector{Vector{UInt32}} = []

    for i in 1:number_encryption_runs
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, i)
        if i == 1
            add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, i)
        end
        add_starting_constant_values(variables, factors, bits_per_cluster, i)
        push!(expected_states, encrypt(key, nonce, counter))
        counter = counter + 0x1
    end

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()
    # println(number_encryption_runs)
    add_equality_between_keys_and_nonces(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors)
    add_adds_between_counters(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors, additional_variables)

    variables_by_round[1] = union(variables_by_round[1], additional_variables)
    factors_by_round[1] = union(factors_by_round[1], additional_factors)

    # Need to pass messages through the factor graph to reach the correct values
    belief_propogate_through_graph_forwards(variables, factors, variables_by_round, factors_by_round, adds_by_round, bits_per_cluster, number_encryption_runs)

    actual_most_likely_states = [[read_most_likely_value_from_variable(variables, string(i, "_", location_execution_counts[i]), bits_per_cluster, encryption_run) for i in 1:16] for encryption_run in 1:number_encryption_runs]
    println(actual_most_likely_states)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_states == actual_most_likely_states
end

function test_vec_3_full_encryption()
    key = zeros(UInt32, 8)
    key[8] = 0x01000000
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 1

    bits_per_cluster = 2
    number_encryption_runs = 3
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    expected_states::Vector{Vector{UInt32}} = []

    for i in 1:number_encryption_runs
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, i)
        if i == 1
            add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, i)
        end
        add_starting_constant_values(variables, factors, bits_per_cluster, i)
        push!(expected_states, encrypt(key, nonce, counter))
        counter = counter + 0x1
    end

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()
    # println(number_encryption_runs)
    add_equality_between_keys_and_nonces(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors)
    add_adds_between_counters(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors, additional_variables)

    variables_by_round[1] = union(variables_by_round[1], additional_variables)
    factors_by_round[1] = union(factors_by_round[1], additional_factors)

    # Need to pass messages through the factor graph to reach the correct values
    belief_propogate_through_graph_forwards(variables, factors, variables_by_round, factors_by_round, adds_by_round, bits_per_cluster, number_encryption_runs)

    actual_most_likely_states = [[read_most_likely_value_from_variable(variables, string(i, "_", location_execution_counts[i]), bits_per_cluster, encryption_run) for i in 1:16] for encryption_run in 1:number_encryption_runs]
    println(actual_most_likely_states)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_states == actual_most_likely_states
end

function test_vec_4_full_encryption()
    key = zeros(UInt32, 8)
    key[1] = 0x0000ff00
    nonce = zeros(UInt32, 3)
    counter::UInt32 = 2

    bits_per_cluster = 2
    number_encryption_runs = 2
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    expected_states::Vector{Vector{UInt32}} = []

    for i in 1:number_encryption_runs
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, i)
        if i == 1
            add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, i)
        end
        add_starting_constant_values(variables, factors, bits_per_cluster, i)
        push!(expected_states, encrypt(key, nonce, counter))
        counter = counter + 0x1
    end

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()
    # println(number_encryption_runs)
    add_equality_between_keys_and_nonces(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors)
    add_adds_between_counters(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors, additional_variables)

    variables_by_round[1] = union(variables_by_round[1], additional_variables)
    factors_by_round[1] = union(factors_by_round[1], additional_factors)

    # Need to pass messages through the factor graph to reach the correct values
    belief_propogate_through_graph_forwards(variables, factors, variables_by_round, factors_by_round, adds_by_round, bits_per_cluster, number_encryption_runs)

    actual_most_likely_states = [[read_most_likely_value_from_variable(variables, string(i, "_", location_execution_counts[i]), bits_per_cluster, encryption_run) for i in 1:16] for encryption_run in 1:number_encryption_runs]
    println(actual_most_likely_states)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_states == actual_most_likely_states
end

function test_vec_5_full_encryption()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    nonce[3] = 0x02000000
    counter::UInt32 = 0

    bits_per_cluster = 2
    number_encryption_runs = 2
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    expected_states::Vector{Vector{UInt32}} = []

    for i in 1:number_encryption_runs
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, i)
        if i == 1
            add_distribution_of_initial_values(variables, factors, bits_per_cluster, key, nonce, counter, i)
        end
        add_starting_constant_values(variables, factors, bits_per_cluster, i)
        push!(expected_states, encrypt(key, nonce, counter))
        counter = counter + 0x1
    end

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()
    # println(number_encryption_runs)
    add_equality_between_keys_and_nonces(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors)
    add_adds_between_counters(variables, factors, bits_per_cluster, number_encryption_runs, additional_factors, additional_variables)

    variables_by_round[1] = union(variables_by_round[1], additional_variables)
    factors_by_round[1] = union(factors_by_round[1], additional_factors)

    # Need to pass messages through the factor graph to reach the correct values
    belief_propogate_through_graph_forwards(variables, factors, variables_by_round, factors_by_round, adds_by_round, bits_per_cluster, number_encryption_runs)

    actual_most_likely_states = [[read_most_likely_value_from_variable(variables, string(i, "_", location_execution_counts[i]), bits_per_cluster, encryption_run) for i in 1:16] for encryption_run in 1:number_encryption_runs]
    println(actual_most_likely_states)
    # Check that the most likely values of the output are equal to those of the input
    @test expected_states == actual_most_likely_states
end

# I think that a good strategy for ordering the passing of information around would be 
# to have a priority queue, where the variables with the highest change in entropy are first,
# and then select the highest and do the message passing for it and its neighbours recomputing
# the changes in entropy accordingly updating the priority queue allowing for a dynamic moving of
# information around the graph (may potentially want some additional structure for triggering a better
# passing through the and trees but may not be required.  Then can take the total entropy of the graph
# every few hundred/thousand updates to see how it is fairing (which could now just be a sum rather than)
# needing to calculate it on every iteration

test_vec_sec_2_full_encryption()
test_vec_1_full_encryption()
test_vec_2_full_encryption()
test_vec_3_full_encryption()
test_vec_4_full_encryption()
test_vec_5_full_encryption()