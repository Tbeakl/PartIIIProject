using Plots

include("node.jl")
include("messages.jl")
include("dynamic_message_scheduling.jl")
include("add_leakage_to_graph.jl")
include("leakage_functions.jl")
include("chacha.jl")
include("chacha_factor_graph.jl")
include("byte_hamming_weight_traces.jl")

# key = zeros(UInt32, 8)
# nonce = zeros(UInt32, 3)
# counter::UInt32 = 0

key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
nonce = [0x09000000, 0x4a000000, 0x00000000]
counter::UInt32 = 1

encryption_trace = encrypt_collect_trace(key, nonce, counter, byte_hamming_weight_for_value)
encryption_output = encrypt(key, nonce, counter)
number_of_bits = 2
hamming_position_table = table_for_hamming_values(number_of_bits)
add_byte_hamming_weight_to_variable = byte_hamming_weight_value_to_function(hamming_position_table )
variables = Dict{String, Variable}()
factors = Dict{String, Factor}()
variables_by_round::Vector{Set{String}} = []
factors_by_round::Vector{Set{String}} = []
adds_by_round::Vector{Vector{Int64}} = []
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
add_starting_constant_values(variables, factors, number_of_bits)
add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter)

add_initial_key_dist(variables, factors, number_of_bits, byte_hamming_weight_for_value.(key), add_byte_hamming_weight_to_variable)
# Need to add some noisy distributions to the initial values for the counters to see how that does

for i in 1:16
    set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits)
end
println("Starting to add the factors for the trace")
add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, add_byte_hamming_weight_to_variable)
println("Added the factors for the trace")

# # Perform the dynamic message passing around the graph to push information through the entire graph
# dynamic_belief_propogate_through_graph(variables, factors, 300_000)

# for i in 1:8
#     println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits), base=16))
# end

# This seems to heavily converge to bad conditions, may be better to do some more full passes at some point to hopefully get to a correct convergence

internal_variables = Set{String}()
internal_factors = Set{String}()
for (i, vars_for_round) in enumerate(variables_by_round)
    union!(internal_variables, vars_for_round)
    union!(internal_factors, factors_by_round[i])
end
all_adds = 0:adds_by_round[end][end]
# I think there is a mistake in the way that the key is put into the factor graph because with reduced noise levels it all goes to zeros
# therefore will look at implementing a simple none noisy version of the calculations and adding of information because that should definitely
# not go too all zeros, because it could be some problem with the way the noise is incorporated into the calculations

# Do a few passes of full updates across everything to attempt to get to better convergence of the result

for (fact_name, factor) in factors
    factor_to_variable_messages(factor)
end

entropy_in_graph::Vector{Float64} = Vector{Float64}()
# With 4 bits per cluster approximately 9/10 per minute (after optimisations about 17 rounds a minute)
# With 2 bits per cluster approximately 28.6 iterations per minute (after some optimisations around 75 a minute maybe slightly more because time taken to make graph included in that)
for i in 1:80
    println(i)
    # I think it is slightly better to go forwards and backwards through the graph
    # compared to just doing at random but it does not seem to be a big help
    # and it does not appear that helping to speed up going through adds is massicvely useful
    belief_propagate_forwards_and_back_through_graph(variables, factors, variables_by_round, factors_by_round, 1)
    
    push!(entropy_in_graph, total_entropy_of_graph(variables))
    println(entropy_in_graph[end])
end

plot(entropy_in_graph)
#31425.174177981593
#31572.453308735967