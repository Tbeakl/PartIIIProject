using Plots, Base.Threads, Random, NPZ
include("../../belief_propagation/node.jl")
gr()
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("template_attack_traces.jl")

number_of_bits::Int64 = 2
number_of_encryption_traces::Int64 = 1

dimensions::Int64 = 8
signal_to_noise_ratio::Float64 = 1.4
key::Vector{UInt32} = generate_random_key()
nonce::Vector{UInt32} = generate_random_nonce()
counter::UInt32 = generate_random_counter()

# The factor graph has clearly been broken because the adds are not working around the graph

noise = noise_distribution_fixed_standard_dev(1.0, dimensions)
mean_vectors = generate_mean_vectors(noise, signal_to_noise_ratio, 255)
add_byte_key_template_to_variable = byte_template_value_to_function(mean_vectors, noise)
add_byte_intermediate_add_template_to_variable = byte_template_value_to_function(mean_vectors, noise)
add_byte_intermediate_rot_template_to_variable = byte_template_value_to_function(mean_vectors, noise)

# base_path_key_mean_vectors =              "D:\\ChaChaData\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_A00\\template_expect_b"
# base_path_intermediate_add_mean_vectors = "D:\\ChaChaData\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_B00\\template_expect_b"
# base_path_intermediate_rot_mean_vectors = "D:\\ChaChaData\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_B00\\template_expect_b"
# add_byte_key_template_to_variable = byte_template_path_to_function(base_path_key_mean_vectors, noise, 200)
# add_byte_intermediate_add_template_to_variable = byte_template_path_to_function(base_path_intermediate_add_mean_vectors, noise, 40)
# add_byte_intermediate_rot_template_to_variable = byte_template_path_to_function(base_path_intermediate_rot_mean_vectors, noise, 40)

variables = Dict{String,AbsVariable}()
factors = Dict{String,AbsFactor}()
variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
for encryption_run_number in 1:number_of_encryption_traces
    encryption_trace = encrypt_collect_trace(key, nonce, counter, byte_values_for_input)
    encryption_output = encrypt(key, nonce, counter)
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
    if encryption_run_number == 1
        add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
    end
    # add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter, 1)
    add_initial_nonce_and_counter_dist(variables, factors, number_of_bits, byte_values_for_input.(nonce), byte_values_for_input(counter), encryption_run_number, add_byte_key_template_to_variable)
    add_initial_key_dist(variables, factors, number_of_bits, byte_values_for_input.(key), encryption_run_number, add_byte_key_template_to_variable)
    # Need to add some noisy distributions to the initial values for the counters to see how that does

    for i in 1:16
        add_byte_key_template_to_variable(byte_values_for_input(encryption_output[i]), variables, factors, number_of_bits, string(i, "_", location_execution_counts[i]), encryption_run_number)
        # set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits, 1)
    end
    println("Starting to add the factors for the trace")
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, encryption_run_number, add_byte_intermediate_add_template_to_variable, add_byte_intermediate_rot_template_to_variable)
    println("Added the factors for the trace")
    counter = counter + 0x1
end

additional_variables::Set{String} = Set{String}()
additional_factors::Set{String} = Set{String}()

if number_of_encryption_traces > 1
    # Add constraint between the counters
    add_adds_between_counters(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors, additional_variables)
end

all_variables::Vector{String} = [keys(variables)...]
all_factors::Vector{String} = [keys(factors)...]

heatmap_plotting_function = plot_current_entropy(variables)
visualisation_of_entropy::Vector{Matrix{Float64}} = []
visualisation_variables, x_labels, y_labels = make_positions_to_var_names(number_of_bits, 1)
tot_entropy_over_time::Vector{Float64} = []

# # Perform the dynamic message passing around the graph to push information through the entire graph
# dynamic_belief_propogate_through_graph(variables, factors, 100_000)

# This seems to heavily converge to bad conditions, may be better to do some more full passes at some point to hopefully get to a correct convergence

# There is some issue with the passing of the probabilities in this because they are all going to zeros
# which is a major problem, not really sure where it is coming from because of the need 

Threads.@threads for fact_name in all_factors
    factor_to_variable_messages(factors[fact_name])
end

Threads.@threads for var_name in all_variables
    variable_to_factor_messages(variables[var_name])
end

# Think it could be helpful to initially push through the values which have actually been set away from uniform
# by the leakage, because there are parts which are left unknown which may not be ideal to push through immediately because
# they are providing no information about what the key actually is - I will currently just do this for adds because
# I think these will be the main source of the issues although there could be ones also for the rotations probably should 
# just add the correct leakages in for before the rotation has taken place as well because that might actually be recoverable

# for add_nums in adds_by_round
#     for add_num in add_nums
#         for encryption_run in 1:number_of_encryption_traces
#             belief_propagate_through_add(variables, factors, number_of_bits, add_num, encryption_run)
#         end
#     end
# end

update_all_entropies(variables, all_variables)
push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))
push!(tot_entropy_over_time, total_entropy_of_graph(variables))
println(tot_entropy_over_time[end])

internal_factors = [union(additional_factors, factors_by_round[:]...)...]
internal_variables = [union(additional_variables, variables_by_round[:]...)...]

initial_number_of_iterations = 200

for i in 1:initial_number_of_iterations
    println(i)
    # timing_info = @timed begin
        Threads.@threads for var_name in internal_variables
            variable_to_factor_messages(variables[var_name], 0.8)
        end
        Threads.@threads for fact_name in internal_factors
            factor_to_variable_messages(factors[fact_name], 0.8)
        end
    # end
    # push!(execution_times, timing_info.time)
    # push!(gc_execution_times, timing_info.gctime)
    update_all_entropies(variables, all_variables)
    push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])
    if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6
        break
    end
end

# println("Average time: ", mean(execution_times), " Standard deviation: ", std(execution_times))
# println("Average gc time: ", mean(gc_execution_times), " Standard deviation: ", std(gc_execution_times))

# Potentially could look at just the very start and end for checking for if we have reached convergence because we 
# don't really need to solve the centre of the graph if the start and end are correct and will not really change any more
# number_of_end_rounds_to_check = 2
# vars_to_check = [union(variables_by_round[begin:number_of_end_rounds_to_check]..., variables_by_round[end-number_of_end_rounds_to_check-1:end]...)...]
# calculate_entropy_of_ends() = sum([variables[var].current_entropy for var in vars_to_check])


# number_of_iterations_of_ends = 200 # 200
# rounds_for_ends = 5
# variables_at_ends = [union(additional_variables, variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
# factors_at_ends = [union(additional_factors, factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]
# prev_ent = total_entropy_of_graph(variables)
# for i in 1:number_of_iterations_of_ends
#     println(i)
#     Threads.@threads for var_name in variables_at_ends
#         variable_to_factor_messages(variables[var_name])
#     end
#     Threads.@threads for fact_name in factors_at_ends
#         factor_to_variable_messages(factors[fact_name])
#     end
#     update_all_entropies(variables, variables_at_ends)
#     push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

#     push!(tot_entropy_over_time, total_entropy_of_graph(variables))
#     println(tot_entropy_over_time[end])
#     cur_ent = total_entropy_of_graph(variables)
#     if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6 || abs(cur_ent - prev_ent) <= 1e-7
#         break
#     end
#     prev_ent = cur_ent
# end

read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1) for i in 1:8]
for i in 1:8
    println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1), base=16))
end

plot(tot_entropy_over_time)

anim = @animate for i in eachindex(visualisation_of_entropy)
    heatmap(visualisation_of_entropy[i]; title=string("Round ", i - 1, " entropy of variables"), clim=(0, number_of_bits)) # 
end
# heatmap(visualisation_of_entropy[1]; title=string("Round ", 0, " entropy of variables")) # clim=(0, number_of_bits),
gif(anim, "test.gif", fps=20)