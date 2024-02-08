using Plots, Base.Threads, Random
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("template_attack_traces.jl")
Random.seed!(1234)
# None converging values
# 0xd13796db
# 0xe1223fde
# 0x1d4d42b5
# 0x95b5b6da
# 0x7d62bb82
# 0x7e15fa8c
# 0xe65c0736
# 0xb0bbe3bb

# 0xbaab5977
# 0x0071efdd
# 0xf2863e8e

# 0x9d6f2457

# key = zeros(UInt32, 8)
# nonce = zeros(UInt32, 3)
# counter::UInt32 = 0

# key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
# nonce = [0x09000000, 0x4a000000, 0x00000000]
# counter::UInt32 = 1


number_of_bits = 2
generate_fresh_graph = true
file_to_load_or_store = "base_graph.blob"

dimensions = 4
signal_to_noise_ratio = 1.7
key = generate_random_key()
nonce = generate_random_nonce()
counter = generate_random_counter()

encryption_trace = encrypt_collect_trace(key, nonce, counter, byte_values_for_input)
encryption_output = encrypt(key, nonce, counter)

# noise = noise_distribution_random_variances(dimensions)
noise = noise_distribution_fixed_standard_dev(0.4)
distribution_from_hamming_weight_per_bit = noise_distribution_fixed_standard_dev(0.1)
# mean_vectors = generate_mean_vectors(noise, signal_to_noise_ratio, 255)
mean_vectors = generate_mean_vectors_based_on_hamming_weights(distribution_from_hamming_weight_per_bit, 8)
add_byte_template_to_variable = byte_template_value_to_function(mean_vectors, noise)

variables = Dict{String,Variable{Factor}}()
factors = Dict{String,Factor{Variable}}()
variables_by_round::Vector{Set{String}} = []
factors_by_round::Vector{Set{String}} = []
adds_by_round::Vector{Vector{Int64}} = []
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
add_starting_constant_values(variables, factors, number_of_bits)
add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter)
add_initial_key_dist(variables, factors, number_of_bits, byte_values_for_input.(key), add_byte_template_to_variable)
# Need to add some noisy distributions to the initial values for the counters to see how that does

for i in 1:16
    set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits)
end
println("Starting to add the factors for the trace")
add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, add_byte_template_to_variable)
println("Added the factors for the trace")

all_variables = [keys(variables)...]
all_factors = [keys(factors)...]

heatmap_plotting_function = plot_current_entropy(variables)
visualisation_of_entropy::Vector{Matrix{Float64}} = []
visualisation_variables, x_labels, y_labels = make_positions_to_var_names(number_of_bits)
tot_entropy_over_time::Vector{Float64} = []

# # Perform the dynamic message passing around the graph to push information through the entire graph
# dynamic_belief_propogate_through_graph(variables, factors, 100_000)

# This seems to heavily converge to bad conditions, may be better to do some more full passes at some point to hopefully get to a correct convergence

# There is some issue with the passing of the probabilities in this because they are all going to zeros
# which is a major problem, not really sure where it is coming from because of the need 

Threads.@threads for fact_name in all_factors
    factor_to_variable_messages(factors[fact_name])
end

# Think it could be helpful to initially push through the values which have actually been set away from uniform
# by the leakage, because there are parts which are left unknown which may not be ideal to push through immediately because
# they are providing no information about what the key actually is - I will currently just do this for adds because
# I think these will be the main source of the issues although there could be ones also for the rotations probably should 
# just add the correct leakages in for before the rotation has taken place as well because that might actually be recoverable

for i in 0:adds_by_round[end][end]
    belief_propagate_through_add(variables, factors, number_of_bits, i)
end

update_all_entropies(variables, all_variables)
push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))
push!(tot_entropy_over_time, total_entropy_of_graph(variables))
println(tot_entropy_over_time[end])

internal_factors = [union(factors_by_round[:]...)...]
internal_variables = [union(variables_by_round[:]...)...]

initial_number_of_iterations = 20

for i in 1:initial_number_of_iterations
    println(i)
    Threads.@threads for var_name in internal_variables
        variable_to_factor_messages(variables[var_name])
    end
    Threads.@threads for fact_name in internal_factors
        factor_to_variable_messages(factors[fact_name])
    end
    update_all_entropies(variables, all_variables)
    push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))
    
    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])
    if tot_entropy_over_time[end] < 1e-3 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-3
        break
    end
end

# Potentially could look at just the very start and end for checking for if we have reached convergence because we 
# don't really need to solve the centre of the graph if the start and end are correct and will not really change any more
number_of_end_rounds_to_check = 2
vars_to_check = [union(variables_by_round[begin:number_of_end_rounds_to_check]..., variables_by_round[end- number_of_end_rounds_to_check - 1:end]...)...]
calculate_entropy_of_ends() = sum([variables[var].current_entropy for var in vars_to_check])


number_of_iterations_of_ends = 400
rounds_for_ends = 5
variables_at_ends = [union(variables_by_round[begin:rounds_for_ends]..., variables_by_round[end- rounds_for_ends - 1:end]...)...]
factors_at_ends = [union(factors_by_round[begin:rounds_for_ends]..., factors_by_round[end- rounds_for_ends - 1:end]...)...]
prev_ent = calculate_entropy_of_ends()
for i in 1:number_of_iterations_of_ends
    println(i)
    Threads.@threads for var_name in variables_at_ends
        variable_to_factor_messages(variables[var_name])
    end
    Threads.@threads for fact_name in factors_at_ends
        factor_to_variable_messages(factors[fact_name])
    end
    update_all_entropies(variables, variables_at_ends)
    push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])
    cur_ent = calculate_entropy_of_ends()
    if tot_entropy_over_time[end] < 1e-3 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-3 || abs(cur_ent - prev_ent) <= 1e-7
        break
    end
    prev_ent = cur_ent
end

read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits) for i in 1:8]
for i in 1:8
    println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits), base=16))
end

plot(tot_entropy_over_time)

anim = @animate for i in 1:length(visualisation_of_entropy)
    heatmap(visualisation_of_entropy[i]; title=string("Round ", i - 1, " entropy of variables"), clim=(0, number_of_bits)) # 
end
# heatmap(visualisation_of_entropy[1]; title=string("Round ", 0, " entropy of variables")) # clim=(0, number_of_bits),
gif(anim, "0_85_end_rounds_damping_1234.gif", fps=30)