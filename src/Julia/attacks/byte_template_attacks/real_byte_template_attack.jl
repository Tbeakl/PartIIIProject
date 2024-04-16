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
number_of_encryption_traces::Int64 = 1
number_of_values_averaged_over_key_leakage = 40

base_path_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/initial_templates/"
add_byte_template_function = real_byte_template_path_to_function(base_path_templates)

base_key_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/initial_templates/"

variables = Dict{String,AbsVariable}()
factors = Dict{String,AbsFactor}()
variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]

key_number = 10
key, nonce, counter, encryption_trace = load_attack_trace(key_number, 1)
for encryption_run_number in 1:number_of_encryption_traces
    key, nonce, counter, encryption_trace = load_attack_trace(key_number, encryption_run_number - 1)
    encryption_output = encrypt(key, nonce, counter)
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
    add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)

    # Here we are assuming that we have known nonce and counter
    add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter, 1)

    println("Starting to add the factors for the trace")
    add_leakage_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, encryption_run_number, add_byte_template_function)
    println("Added the factors for the trace")
    println("Starting adding key distribution")
    # add_initial_key_distribution_from_leakage_trace(encryption_trace, variables, factors, number_of_bits, encryption_run_number, base_key_templates)
    add_initial_key_distribution_from_simulated_leakage(byte_values_for_input.(key), variables, factors, number_of_bits, encryption_run_number, base_key_templates, number_of_values_averaged_over_key_leakage)
    println("Added key distribution")

    # Also want to add in the known outputs of the encryption as part of the model
    for i in 1:16
        set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits, encryption_run_number)
    end
end

additional_variables::Set{String} = Set{String}()
additional_factors::Set{String} = Set{String}()

if number_of_encryption_traces > 1
    # Add equality constraint between all the encryption runs and the add constraint between the counters
    add_equality_between_keys_and_nonces(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors)
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

for add_nums in adds_by_round
    for add_num in add_nums
        for encryption_run in 1:number_of_encryption_traces
            belief_propagate_through_add(variables, factors, number_of_bits, add_num, encryption_run)
        end
    end
end

update_all_entropies(variables, all_variables)
push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))
push!(tot_entropy_over_time, total_entropy_of_graph(variables))
println(tot_entropy_over_time[end])

internal_factors = [union(additional_factors, factors_by_round[:]...)...]
internal_variables = [union(additional_variables, variables_by_round[:]...)...]

initial_number_of_iterations = 100

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
    if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6
        break
    end
end

# number_of_iterations_of_ends = 200
# rounds_for_ends = 2
# variables_at_ends = [union(additional_variables, variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
# factors_at_ends = [union(additional_factors, factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]
# for i in 1:number_of_iterations_of_ends
#     println(i)
#     Threads.@threads for var_name in variables_at_ends
#         variable_to_factor_messages(variables[var_name], .8)
#     end
#     Threads.@threads for fact_name in factors_at_ends
#         factor_to_variable_messages(factors[fact_name], .8)
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

plot(tot_entropy_over_time)

anim = @animate for i in eachindex(visualisation_of_entropy)
    heatmap(visualisation_of_entropy[i]; title=string("Round ", i - 1, " entropy of variables"), clim=(0, number_of_bits)) # 
end
# heatmap(visualisation_of_entropy[1]; title=string("Round ", 0, " entropy of variables")) # clim=(0, number_of_bits),
gif(anim, "test.gif", fps=10)