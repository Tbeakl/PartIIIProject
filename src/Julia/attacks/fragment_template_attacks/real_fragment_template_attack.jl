using Plots, Base.Threads, Random, HDF5
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("../../encryption/rank_estimation.jl")
include("template_attack_traces.jl")

# initial_key_number = parse(Int64, ARGS[1])
# final_key_number = parse(Int64, ARGS[2])

number_of_bits::Int64 = 8
bits_per_template::Int64 = 16
dimensions_per_template::Int64 = 16
number_of_encryption_traces::Int64 = 15
number_of_values_averaged_over_key_leakage = 1
number_of_values_trace_averaged_over = 1

base_path_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/initial_templates_sixteen_bit_templates/sparse_50_detailed_50/"
base_key_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/initial_templates_sixteen_bit_templates/sparse_50_detailed_50/"
base_trace_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/captures/ChaChaRecordings_2/recording_attack_counter_constant_"

base_evaluation_of_ranks = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/constant_counter_4_16/"
key_number = 1
# for key_number in initial_key_number:final_key_number
    variables = Dict{String,AbsVariable}()
    factors = Dict{String,AbsFactor}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
    println("Key number: ", key_number)
    all_traces = zeros(Float32, number_of_encryption_traces, 14989)
    key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, 1)
    for i in 1:number_of_encryption_traces
        key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, i - 1)
        all_traces[i, :] = encryption_trace
    end
    encryption_output = encrypt(key, nonce, counter)

    for encryption_run_number in 1:1
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
        if encryption_run_number == 1
            add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
        end

        # Here we are assuming that we have known nonce and counter
        add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter, 1)

        println("Starting to add the factors for the trace")
        add_byte_template_function = real_byte_template_path_to_function(base_path_templates, bits_per_template, dimensions_per_template, all_traces)
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
        bits_per_template,
        dimensions_per_template,
        all_traces)
    println("Added key distribution")

    additional_variables::Set{String} = Set{String}()
    additional_factors::Set{String} = Set{String}()

    all_variables::Vector{String} = [keys(variables)...]
    all_factors::Vector{String} = [keys(factors)...]

    heatmap_plotting_function = plot_current_entropy(variables)
    visualisation_of_entropy::Vector{Matrix{Float64}} = []
    visualisation_variables, x_labels, y_labels = make_positions_to_var_names(number_of_bits, 1)
    tot_entropy_over_time::Vector{Float64} = []

    Threads.@threads for fact_name in all_factors
        factor_to_variable_messages(factors[fact_name])
    end

    Threads.@threads for var_name in all_variables
        variable_to_factor_messages(variables[var_name])
    end


    initial_likelihood_tables = make_log_likelihood_tables_for_key(variables, number_of_bits)
    # Replace all the -Inf with like -800 because of the issues associated with expoentiatation
    initial_estimated_rank = rank_estimate_key(key, initial_likelihood_tables, 500, number_of_bits)

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
            variable_to_factor_messages(variables[var_name], 0.99)
        end
        Threads.@threads for fact_name in internal_factors
            factor_to_variable_messages(factors[fact_name], 0.99)
        end
        update_all_entropies(variables, all_variables)
        push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

        push!(tot_entropy_over_time, total_entropy_of_graph(variables))
        println(tot_entropy_over_time[end])
        if tot_entropy_over_time[end] < 1 || tot_entropy_over_time[end-1] - tot_entropy_over_time[end] <= 0.5
            break
        end
    end

    final_likelihood_tables = make_log_likelihood_tables_for_key(variables, number_of_bits)
    # Replace all the -Inf with like -800 because of the issues associated with expoentiatation
    final_estimated_rank = rank_estimate_key(key, final_likelihood_tables, 500, number_of_bits)

    read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1) for i in 1:8]
    for i in 1:8
        println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1), base=16))
    end

    result_fid = h5open(string(base_evaluation_of_ranks, key_number, ".hdf5"), "w")
    result_fid["initial_likelihood_tables"] = reduce(vcat,transpose.(initial_likelihood_tables))
    initial_estimated_rank = max(1, initial_estimated_rank)
    result_fid["initial_estimated_rank_log2"] = Float64(log2(initial_estimated_rank))
    result_fid["final_likelihood_tables"] = reduce(vcat,transpose.(final_likelihood_tables))
    final_estimated_rank = max(1, final_estimated_rank)
    result_fid["final_estimated_rank_log2"] = Float64(log2(final_estimated_rank))
    result_fid["entropy_over_time"] = tot_entropy_over_time
    close(result_fid)
# end

# anim = @animate for i in eachindex(visualisation_of_entropy) #(length(visualisation_of_entropy)-1):length(visualisation_of_entropy)#
#     Plots.heatmap(visualisation_of_entropy[i]; title=string("Round ", i - 1, " entropy of variables"), clim=(0, number_of_bits)) # 
# end
# # heatmap(visualisation_of_entropy[1]; title=string("Round ", 0, " entropy of variables")) # clim=(0, number_of_bits),
# gif(anim, "test.gif", fps=10)