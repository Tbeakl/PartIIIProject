using Plots, Base.Threads, Random, HDF5, Distributed
@everywhere include("../../belief_propagation/node.jl")
@everywhere include("../../belief_propagation/messages.jl")
@everywhere include("../../chacha_factor_graph/chacha_factor_graph.jl")
@everywhere include("../../chacha_factor_graph/add_leakage_to_graph.jl")
@everywhere include("../../chacha_factor_graph/heatmap_visualisation.jl")
@everywhere include("../../encryption/leakage_functions.jl")
@everywhere include("../../encryption/chacha.jl")
@everywhere include("../../encryption/rank_estimation.jl")
@everywhere include("template_attack_traces.jl")

@everywhere function evaluate_key_number(number_of_bits::Int64,
    bits_per_template::Int64,
    dimensions_per_template::Int64,
    number_of_encryption_traces::Int64,
    key_number::Int64)

    path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

    base_path_templates = path_to_data * "attack_profiling/32_volatile/initial_templates_8bits/"
    base_key_templates = path_to_data * "attack_profiling/32_volatile/initial_templates_8bits/"
    base_trace_path = path_to_data * "captures/ChaChaRecordings_3/recording_attack_counter_from_random_"

    base_evaluation_of_ranks = path_to_data * "evaluation/incremented_random_counter_32_volatile_10_8/"

    path_to_mean = path_to_data * "attack_profiling/32_volatile/mean_trace.hdf5"

    damping_factor::Float64 = 0.95

    result_path = string(base_evaluation_of_ranks, key_number, ".hdf5")
    if !ispath(result_path)
        variables = Dict{String,AbsVariable}()
        factors = Dict{String,AbsFactor}()
        variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
        factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
        adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
        println("Key number: ", key_number)
        all_traces = zeros(Float32, number_of_encryption_traces, 129960)
        key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, 1, path_to_mean)
        for i in 1:number_of_encryption_traces
            key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, i - 1, path_to_mean)
            all_traces[i, :] = encryption_trace
        end
        encryption_output = encrypt(key, nonce, counter)

        for encryption_run_number in 1:number_of_encryption_traces
            key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, encryption_run_number - 1, path_to_mean)
            encryption_output = encrypt(key, nonce, counter)
            location_execution_counts = zeros(Int64, 16)
            chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
            if encryption_run_number == 1
                add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
            end

            # Here we are assuming that we have known nonce and counter
            add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter, encryption_run_number)
            # add_initial_counter_distribution_from_leakage_traces_set_of_templates(variables,
            #     factors,
            #     number_of_bits,
            #     base_key_templates,
            #     bits_per_template,
            #     dimensions_per_template,
            #     encryption_trace,
            #     encryption_run_number)

            println("Starting to add the factors for the trace")
            add_byte_template_function = real_byte_template_path_to_function(base_path_templates, bits_per_template, dimensions_per_template, encryption_trace')
            add_leakage_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, encryption_run_number, add_byte_template_function)
            println("Added the factors for the trace")

            # Also want to add in the known outputs of the encryption as part of the model
            for i in 1:16
                set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits, encryption_run_number)
            end
            # all_distribution_of_output(
            #     variables,
            #     factors,
            #     number_of_bits,
            #     base_key_templates,
            #     bits_per_template,
            #     dimensions_per_template,
            #     location_execution_counts,
            #     encryption_trace,
            #     encryption_run_number
            # )
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

        # add_initial_nonce_distribution_from_leakage_traces_set_of_templates(
        #     variables,
        #     factors,
        #     number_of_bits,
        #     base_key_templates,
        #     bits_per_template,
        #     dimensions_per_template,
        #     all_traces)
        println("Added key distribution")

        additional_variables::Set{String} = Set{String}()
        additional_factors::Set{String} = Set{String}()

        all_variables::Vector{String} = [keys(variables)...]
        all_factors::Vector{String} = [keys(factors)...]

        heatmap_plotting_function = plot_current_entropy(variables)
        visualisation_of_entropy::Vector{Matrix{Float64}} = []
        visualisation_variables, x_labels, y_labels = make_positions_to_var_names(number_of_bits, 1)
        tot_entropy_over_time::Vector{Float64} = []
        # println("Here")
        for fact_name in all_factors
            factor_to_variable_messages(factors[fact_name])
        end

        for var_name in all_variables
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

        initial_number_of_iterations = 20

        for i in 1:initial_number_of_iterations
            println(i)
            for var_name in internal_variables #Threads.@threads 
                variable_to_factor_messages(variables[var_name], damping_factor)
            end
            for fact_name in internal_factors
                factor_to_variable_messages(factors[fact_name], damping_factor)
            end
            update_all_entropies(variables, all_variables)
            push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

            push!(tot_entropy_over_time, total_entropy_of_graph(variables))
            println(tot_entropy_over_time[end])
            if tot_entropy_over_time[end] < 100 || tot_entropy_over_time[end-1] - tot_entropy_over_time[end] <= 0.5 || isnan(tot_entropy_over_time[end])
                break
            end
        end

        number_of_iterations_of_ends = 200
        rounds_for_ends = 2
        variables_at_ends = [union(additional_variables, variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
        factors_at_ends = [union(additional_factors, factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]
        for i in 1:number_of_iterations_of_ends
            println(i)
            for var_name in variables_at_ends
                variable_to_factor_messages(variables[var_name], damping_factor)
            end
            for fact_name in factors_at_ends
                factor_to_variable_messages(factors[fact_name], damping_factor)
            end
            update_all_entropies(variables, variables_at_ends)
            push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

            push!(tot_entropy_over_time, total_entropy_of_graph(variables))
            println(tot_entropy_over_time[end])
            if tot_entropy_over_time[end] < 100 || tot_entropy_over_time[end-1] - tot_entropy_over_time[end] <= 0.5 || isnan(tot_entropy_over_time[end])
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

        result_fid = h5open(result_path, "w")
        result_fid["initial_likelihood_tables"] = reduce(vcat, transpose.(initial_likelihood_tables))
        initial_estimated_rank = max(1, initial_estimated_rank)
        result_fid["initial_estimated_rank_log2"] = Float64(log2(initial_estimated_rank))
        result_fid["final_likelihood_tables"] = reduce(vcat, transpose.(final_likelihood_tables))
        final_estimated_rank = max(1, final_estimated_rank)
        result_fid["final_estimated_rank_log2"] = Float64(log2(final_estimated_rank))
        result_fid["entropy_over_time"] = tot_entropy_over_time
        close(result_fid)
        read_off_key == key
    else
        1
    end
end

initial_key_number = parse(Int64, ARGS[1])
final_key_number = parse(Int64, ARGS[2])
number_of_bits::Int64 = 8
bits_per_template::Int64 = 8
dimensions_per_template::Int64 = 8
number_of_encryption_traces::Int64 = 10

# trace_numbers = [2, 6, 10, 16, 21, 30, 35, 45, 46, 48, 49, 51, 55, 58, 59, 65, 74, 75, 76, 85, 86, 87, 93, 99, 100,
#     101, 102, 108, 109, 110, 113, 119, 124, 127, 129, 132, 135, 138, 143, 148, 150, 152, 154, 160, 161,
#     169, 172, 179, 180, 181, 186, 188, 192, 200, 206, 211, 213, 215, 216, 218, 221, 236, 238, 242, 245,
#     248, 250, 256, 258, 262, 263, 267, 272, 279, 281, 292, 296, 301, 312, 315, 319, 322, 323, 328, 329,
#     347, 349, 352, 355, 356, 366, 372, 373, 374, 380, 381, 383, 384, 389, 390, 401, 403, 405, 408, 412,
#     415, 416, 426, 431, 432, 439, 441, 442, 444, 447, 448, 450, 455, 457, 461, 471, 472, 477, 478, 488,
#     489, 490, 496, 504, 505, 511, 515, 518, 520, 543, 544, 545, 547, 558, 559, 565, 566, 568, 573, 575,
#     577, 580, 581, 584, 589, 592, 595, 598, 601, 603, 607, 616, 620, 625, 637, 638, 639, 645, 653, 654,
#     659, 662, 663, 665, 674, 677, 682, 689, 690, 692, 700, 701, 705, 710, 716, 721, 722, 723, 728, 735,
#     744, 745, 747, 749, 755, 761, 763, 764, 769, 772, 778, 780, 781, 783, 786, 787, 788, 800, 801, 803,
#     805, 806, 810, 813, 814, 815, 816, 817, 820, 826, 828, 831, 832, 836, 837, 839, 863, 865, 874, 878,
#     880, 883, 884, 890, 891, 892, 893, 897, 899, 900, 902, 911, 913, 918, 919, 920, 923, 929, 932, 934,
#     944, 946, 948, 952, 958, 962, 964, 973, 976, 978, 980, 987, 988, 990, 995, 996, 998]

num_successes = @distributed (+) for i in initial_key_number:final_key_number
    evaluate_key_number(number_of_bits, bits_per_template, dimensions_per_template, number_of_encryption_traces, i)
end
evaluate_key_number(number_of_bits, bits_per_template, dimensions_per_template, number_of_encryption_traces, 2)
println(num_successes)