using Test

include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../belief_propagation/dynamic_message_scheduling.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("exact_value_traces.jl")

function all_zeros_single_bit_clusters()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 1
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)
    add_starting_constant_values(variables, factors, number_of_bits, 1)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter, 1)
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, 1)
    all_var_names = collect(keys(variables))
    for execution_count in 1:3
        for (i,j) in factors
            factor_to_variable_messages(j)
        end
        for (i,j) in variables
            variable_to_factor_messages(j)
        end
        update_all_entropies(variables, all_var_names)
        println("Total entropy after ", execution_count, " passes ", total_entropy_of_graph(variables))
    end
    # Make a copy of all of the distributions to check for a fixed point
    variable_marginal_distributions = Dict{String, Vector{Float64}}()
    for (i,j) in variables
        variable_marginal_distributions[i] = marginal(j)
    end

    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end

    for (i,j) in variables
        @test marginal(j) == variable_marginal_distributions[i]
    end
end

function all_zeros_2_bit_clusters()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 2
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)
    add_starting_constant_values(variables, factors, number_of_bits, 1)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter, 1)
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, 1)
    all_var_names = collect(keys(variables))
    for execution_count in 1:3
        for (i,j) in factors
            factor_to_variable_messages(j)
        end
        for (i,j) in variables
            variable_to_factor_messages(j)
        end
        update_all_entropies(variables, all_var_names)
        println("Total entropy after ", execution_count, " passes ", total_entropy_of_graph(variables))
    end
    # Make a copy of all of the distributions to check for a fixed point
    variable_marginal_distributions = Dict{String, Vector{Float64}}()
    for (i,j) in variables
        variable_marginal_distributions[i] = marginal(j)
    end

    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end

    for (i,j) in variables
        @test marginal(j) == variable_marginal_distributions[i]
    end
end

function all_zeros_4_bit_clusters()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 4
    variables = Dict{String, AbsVariable}()
    factors = Dict{String, AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)
    add_starting_constant_values(variables, factors, number_of_bits, 1)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter, 1)
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, 1)
    all_var_names = collect(keys(variables))
    for execution_count in 1:3
        for (i,j) in factors
            factor_to_variable_messages(j)
        end
        for (i,j) in variables
            variable_to_factor_messages(j)
        end
        update_all_entropies(variables, all_var_names)
        println("Total entropy after ", execution_count, " passes ", total_entropy_of_graph(variables))
    end
    # Make a copy of all of the distributions to check for a fixed point
    variable_marginal_distributions = Dict{String, Vector{Float64}}()
    for (i,j) in variables
        variable_marginal_distributions[i] = marginal(j)
    end

    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end

    for (i,j) in variables
        @test marginal(j) == variable_marginal_distributions[i]
    end
end

all_zeros_single_bit_clusters()
all_zeros_2_bit_clusters()
all_zeros_4_bit_clusters()