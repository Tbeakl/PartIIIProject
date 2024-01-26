using Test
include("chacha_factor_graph.jl")
include("exact_value_traces.jl")
include("dynamic_message_scheduling.jl")

function all_zeros_single_bit_clusters_perturbed_output()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 1
    variables = Dict{String, Variable}()
    factors = Dict{String, Factor}()
    variables_by_round::Vector{Set{String}} = []
    factors_by_round::Vector{Set{String}} = []    
    adds_by_round::Vector{Vector{Int64}} = []
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
    add_starting_constant_values(variables, factors, number_of_bits)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter)
    println(encryption_trace[end])
    # Here flip the lowest bit in the final trace value, which should then allow us to see
    # that information can correctly flow through the graph
    encryption_trace[end] = 0x1 ⊻ encryption_trace[end]
    println(encryption_trace[end])
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits)

    # Perform the dynamic message passing around the graph to push information through the entire graph
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)
    
    # Then need to check that all of the marginal distributions in the graph are all zero
    # because it is an invalid trace
    for (var_name, variable) in variables
        cur_marginal = marginal(variable)
        @test cur_marginal == zeros(length(cur_marginal))
    end
end

function all_zeros_2_bit_clusters_perturbed_output()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 2
    variables = Dict{String, Variable}()
    factors = Dict{String, Factor}()
    variables_by_round::Vector{Set{String}} = []
    factors_by_round::Vector{Set{String}} = []    
    adds_by_round::Vector{Vector{Int64}} = []
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
    add_starting_constant_values(variables, factors, number_of_bits)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter)
    println(encryption_trace[end])
    # Here flip the lowest bit in the final trace value, which should then allow us to see
    # that information can correctly flow through the graph
    encryption_trace[end] = 0x1 ⊻ encryption_trace[end]
    println(encryption_trace[end])
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits)

    # Perform the dynamic message passing around the graph to push information through the entire graph
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)
    
    # Then need to check that all of the marginal distributions in the graph are all zero
    # because it is an invalid trace
    for (var_name, variable) in variables
        cur_marginal = marginal(variable)
        @test cur_marginal == zeros(length(cur_marginal))
    end
end

function all_zeros_4_bit_clusters_perturbed_output()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 4
    variables = Dict{String, Variable}()
    factors = Dict{String, Factor}()
    variables_by_round::Vector{Set{String}} = []
    factors_by_round::Vector{Set{String}} = []    
    adds_by_round::Vector{Vector{Int64}} = []
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
    add_starting_constant_values(variables, factors, number_of_bits)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter)
    println(encryption_trace[end])
    # Here flip the lowest bit in the final trace value, which should then allow us to see
    # that information can correctly flow through the graph
    encryption_trace[end] = 0x1 ⊻ encryption_trace[end]
    println(encryption_trace[end])
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits)

    # Perform the dynamic message passing around the graph to push information through the entire graph
    dynamic_belief_propogate_through_graph(variables, factors, 100_000)
    
    # Then need to check that all of the marginal distributions in the graph are all zero
    # because it is an invalid trace
    for (var_name, variable) in variables
        cur_marginal = marginal(variable)
        @test cur_marginal == zeros(length(cur_marginal))
    end
end

all_zeros_single_bit_clusters_perturbed_output()
all_zeros_2_bit_clusters_perturbed_output()
all_zeros_4_bit_clusters_perturbed_output()