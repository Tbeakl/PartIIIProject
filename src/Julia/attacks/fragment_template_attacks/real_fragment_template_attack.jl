using Plots, Base.Threads, Random, HDF5
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("template_attack_traces.jl")

number_of_bits::Int64 = 8
number_of_encryption_traces::Int64 = 250
number_of_values_averaged_over_key_leakage = 1
number_of_values_trace_averaged_over = 1

base_path_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/initial_templates_HW/sparse_50_detailed_50/"
base_key_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/initial_templates_HW/sparse_50_detailed_50/"
base_trace_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/captures/ChaChaRecordings_2/recording_attack_counter_constant_"

variables = Dict{String,AbsVariable}()
factors = Dict{String,AbsFactor}()
variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]

# testing_fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/profiling_50.hdf5", "r")
# matrix_used = read(testing_fid["downsampled_matrix"])
# intermediate_values_of_profiling = read(testing_fid["intermediate_values"])
# close(testing_fid)

key_number = 1103

all_traces = zeros(Float32, number_of_encryption_traces, 14989)
key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, 1)
for i in 1:number_of_encryption_traces
    key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, i - 1)
    all_traces[i, :] = encryption_trace
end
encryption_output = encrypt(key, nonce, counter)
# There is clearly something not lining up about the intermediate values corresponding to the key or I have really
# badly classified the noise or something, not sure about the other intermediate values because these keys seem way off
elements_of_trace_to_select = append!(repeat([true, false, false, false, true], 320), ones(Bool, 16))
plaintext = zeros(UInt32, 16)
ciphertext = encrypt(key, nonce, counter) .⊻ plaintext
trace = encrypt_collect_trace(key, nonce, counter, full_value)[elements_of_trace_to_select]

all_trace_values = append!(
    key,
    nonce,
    [counter],
    plaintext,
    trace,
    ciphertext)
all_values = UInt32.(collect(Iterators.flatten(all_trace_values)))


intermediate_index_to_plot = 2
template_number_to_use = 4
# profiling_intermediate_value_vector = reduce(vcat, transpose.(byte_values_for_input.(intermediate_values_of_profiling[:, intermediate_index_to_plot])))[:, template_number_to_use]
# traces_to_use_profiling = matrix_used[profiling_intermediate_value_vector.==byte_values_for_input(all_values[intermediate_index_to_plot])[template_number_to_use], :]

p = plot_distribution_of_values_and_means(base_key_templates,
    8, 8,
    intermediate_index_to_plot,
    template_number_to_use, byte_hamming_weight_for_value(all_values[intermediate_index_to_plot])[template_number_to_use],
    all_traces)

for encryption_run_number in 1:1
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
    if encryption_run_number == 1
        add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
    end

    # Here we are assuming that we have known nonce and counter
    add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter, 1)

    println("Starting to add the factors for the trace")
    add_byte_template_function = real_byte_template_path_to_function(base_path_templates, 8, 8, all_traces)
    add_leakage_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, encryption_run_number, add_byte_template_function)
    println("Added the factors for the trace")

    # Also want to add in the known outputs of the encryption as part of the model
    for i in 1:16
        set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits, encryption_run_number)
    end
end

println("Starting adding key distribution")
add_initial_key_distribution_from_leakage_traces_set_of_templates(
    variables,
    factors,
    number_of_bits,
    base_key_templates,
    8,
    8,
    all_traces)
println("Added key distribution")

additional_variables::Set{String} = Set{String}()
additional_factors::Set{String} = Set{String}()

# if number_of_encryption_traces > 1
#     # Add the add constraint between the counters
#     add_adds_between_counters(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors, additional_variables)
# end

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

initial_number_of_iterations = 50

for i in 1:initial_number_of_iterations
    println(i)
    Threads.@threads for var_name in internal_variables
        variable_to_factor_messages(variables[var_name], 0.8)
    end
    Threads.@threads for fact_name in internal_factors
        factor_to_variable_messages(factors[fact_name], 0.8)
    end
    update_all_entropies(variables, all_variables)
    push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])
    if tot_entropy_over_time[end] < 1 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-5
        break
    end
end

# initial_number_of_iterations = 50

# for i in 1:initial_number_of_iterations
#     println(i)
#     Threads.@threads for var_name in internal_variables
#         variable_to_factor_messages(variables[var_name], 0.6)
#     end
#     Threads.@threads for fact_name in internal_factors
#         factor_to_variable_messages(factors[fact_name], 0.6)
#     end
#     update_all_entropies(variables, all_variables)
#     push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

#     push!(tot_entropy_over_time, total_entropy_of_graph(variables))
#     println(tot_entropy_over_time[end])
#     if tot_entropy_over_time[end] < 1 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-5
#         break
#     end
# end


# number_of_iterations_of_ends = 100
# rounds_for_ends = 2
# variables_at_ends = [union(additional_variables, variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
# factors_at_ends = [union(additional_factors, factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]
# # variables_at_ends = [union(variables_by_round[begin+rounds_for_ends:end-rounds_for_ends]...)...]
# # factors_at_ends = [union(factors_by_round[begin+rounds_for_ends:end-rounds_for_ends]...)...]
# for i in 1:number_of_iterations_of_ends
#     println(i)
#     Threads.@threads for var_name in variables_at_ends
#         variable_to_factor_messages(variables[var_name], 0.7)
#     end
#     Threads.@threads for fact_name in factors_at_ends
#         factor_to_variable_messages(factors[fact_name], 0.7)
#     end
#     update_all_entropies(variables, variables_at_ends)
#     push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

#     push!(tot_entropy_over_time, total_entropy_of_graph(variables))
#     println(tot_entropy_over_time[end])
#     if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6
#         break
#     end
# end


# number_of_iterations_of_ends = 200
# rounds_for_ends = 2
# variables_at_ends = [union(additional_variables, variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
# factors_at_ends = [union(additional_factors, factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]
# for i in 1:number_of_iterations_of_ends
#     println(i)
#     Threads.@threads for var_name in variables_at_ends
#         variable_to_factor_messages(variables[var_name], .7)
#     end
#     Threads.@threads for fact_name in factors_at_ends
#         factor_to_variable_messages(factors[fact_name], .7)
#     end
#     update_all_entropies(variables, variables_at_ends)
#     push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

#     push!(tot_entropy_over_time, total_entropy_of_graph(variables))
#     println(tot_entropy_over_time[end])
#     if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6
#         break
#     end
# end


read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1) for i in 1:8]
for i in 1:8
    println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1), base=16))
end

Plots.plot(tot_entropy_over_time)

anim = @animate for i in eachindex(visualisation_of_entropy) #(length(visualisation_of_entropy)-1):length(visualisation_of_entropy)#
    Plots.heatmap(visualisation_of_entropy[i]; title=string("Round ", i - 1, " entropy of variables"), clim=(0, number_of_bits)) # 
end
# heatmap(visualisation_of_entropy[1]; title=string("Round ", 0, " entropy of variables")) # clim=(0, number_of_bits),
gif(anim, "test.gif", fps=10)