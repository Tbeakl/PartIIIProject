include("node.jl")
include("messages.jl")
include("dynamic_message_scheduling.jl")
include("chacha.jl")
include("chacha_factor_graph.jl")
include("noisy_byte_hamming_weight_traces.jl")
include("add_leakage_to_graph.jl")
include("leakage_functions.jl")
# key = zeros(UInt32, 8)
# nonce = zeros(UInt32, 3)
# counter::UInt32 = 0

key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
nonce = [0x09000000, 0x4a000000, 0x00000000]
counter::UInt32 = 1

standard_deviation = 0.7
number_of_bits = 4
hamming_position_table = table_for_hamming_values(number_of_bits)
noise = noise_distribution(standard_deviation)

leakage_function(value) = byte_hamming_weight_for_value(value) .+ rand(noise)

encryption_trace = encrypt_collect_trace(key, nonce, counter, leakage_function)
encryption_output = encrypt(key, nonce, counter)
add_noisy_byte_hamming_weight_to_variable = noisy_byte_hamming_weight_value_to_function(hamming_position_table, noise)

variables = Dict{String, Variable}()
factors = Dict{String, Factor}()
variables_by_round::Vector{Set{String}} = []
factors_by_round::Vector{Set{String}} = []
adds_by_round::Vector{Vector{Int64}} = []
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
add_starting_constant_values(variables, factors, number_of_bits)
add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter)
add_initial_key_dist(variables, factors, number_of_bits, leakage_function.(key), add_noisy_byte_hamming_weight_to_variable)
# Need to add some noisy distributions to the initial values for the counters to see how that does

for i in 1:16
    set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits)
end
println("Starting to add the factors for the trace")
add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, add_noisy_byte_hamming_weight_to_variable)
println("Added the factors for the trace")

# # Perform the dynamic message passing around the graph to push information through the entire graph
# dynamic_belief_propogate_through_graph(variables, factors, 100_000)

# This seems to heavily converge to bad conditions, may be better to do some more full passes at some point to hopefully get to a correct convergence

# There is some issue with the passing of the probabilities in this because they are all going to zeros
# which is a major problem, not really sure where it is coming from because of the need 

for (fact_name, factor) in factors
    factor_to_variable_messages(factor)
end

tot_entropy_over_time::Vector{Float64} = []
push!(tot_entropy_over_time, total_entropy_of_graph(variables))
println(tot_entropy_over_time[end])

for i in 1:100
    println(i)
    belief_propagate_forwards_and_back_through_graph(variables, factors, variables_by_round, factors_by_round, 1)
    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])
    if tot_entropy_over_time[end] < 1e-3
        break
    end
end

for i in 1:8
    println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits), base=16))
end

plot(tot_entropy_over_time)
# I think there is a mistake in the way that the key is put into the factor graph because with reduced noise levels it all goes to zeros
# therefore will look at implementing a simple none noisy version of the calculations and adding of information because that should definitely
# not go too all zeros, because it could be some problem with the way the noise is incorporated into the calculations

# Do a few passes of full updates across everything to attempt to get to better convergence of the result
# for i in 1:30
#     println(i)
#     for fact in internal_factors
#         factor_to_variable_messages(factors[fact])
#     end
#     for var in internal_variables
#         variable_to_factor_messages(variables[var])
#     end
#     println(total_entropy_of_graph(variables))
# end