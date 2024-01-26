include("node.jl")
include("messages.jl")
include("dynamic_message_scheduling.jl")
include("chacha.jl")
include("chacha_factor_graph.jl")
include("noisy_byte_hamming_weight_traces.jl")

key = zeros(UInt32, 8)
nonce = zeros(UInt32, 3)
counter::UInt32 = 0

encryption_trace = encrypt_collect_trace(key, nonce, counter)
encryption_output = encrypt(key, nonce, counter)
number_of_bits = 4
standard_deviation = 0.1
variables = Dict{String, Variable}()
factors = Dict{String, Factor}()
variables_by_round::Vector{Set{String}} = []
factors_by_round::Vector{Set{String}} = []
adds_by_round::Vector{Vector{Int64}} = []
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
add_starting_constant_values(variables, factors, number_of_bits)
add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter)
add_initial_key_dist(variables, factors, number_of_bits, key, standard_deviation)
# Need to add some noisy distributions to the initial values for the counters to see how that does

for i in 1:16
    set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits)
end
println("Starting to add the factors for the trace")
add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, standard_deviation)
println("Added the factors for the trace")

# Perform the dynamic message passing around the graph to push information through the entire graph
dynamic_belief_propogate_through_graph(variables, factors, 800_000)