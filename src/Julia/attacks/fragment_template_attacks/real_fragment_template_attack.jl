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

    base_path_templates = path_to_data * "attack_profiling/8_on_32/initial_templates/"
    base_key_templates = path_to_data * "attack_profiling/8_on_32/initial_templates/"
    base_trace_path = path_to_data * "captures/ChaChaRecordings_8_on_32/recording_attack_counter_from_random_"

    base_evaluation_of_ranks = path_to_data * "evaluation/attack_8_on_32_unknown/"

    damping_factor::Float64 = 0.95

    result_path = string(base_evaluation_of_ranks, key_number, ".hdf5")
    if !ispath(result_path)
        variables = Dict{String,AbsVariable}()
        factors = Dict{String,AbsFactor}()
        variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
        factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
        adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
        println("Key number: ", key_number)
        all_traces = zeros(Float32, number_of_encryption_traces, 399925)
        key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, 1, path_to_data)
        for i in 1:number_of_encryption_traces
            key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, i - 1, path_to_data)
            all_traces[i, :] = encryption_trace
        end
        encryption_output = encrypt(key, nonce, counter)

        for encryption_run_number in 1:number_of_encryption_traces
            location_execution_counts = zeros(Int64, 16)
            chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
            if encryption_run_number == 1
                add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
            end

            # Here we are assuming that we have known nonce and counter
            # add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter, encryption_run_number)
            add_initial_counter_distribution_from_leakage_traces_set_of_templates(variables,
                factors,
                number_of_bits,
                base_key_templates,
                bits_per_template,
                dimensions_per_template,
                encryption_trace,
                encryption_run_number)

            println("Starting to add the factors for the trace")
            add_byte_template_function = real_byte_template_path_to_function(base_path_templates, bits_per_template, dimensions_per_template, all_traces)
            add_leakage_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, encryption_run_number, add_byte_template_function)
            println("Added the factors for the trace")

            # Also want to add in the known outputs of the encryption as part of the model
            # for i in 1:16
            #     set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits, encryption_run_number)
            # end
            all_distribution_of_output(
                variables,
                factors,
                number_of_bits,
                base_key_templates,
                bits_per_template,
                dimensions_per_template,
                location_execution_counts,
                encryption_trace,
                encryption_run_number
            )
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

        add_initial_nonce_distribution_from_leakage_traces_set_of_templates(
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

        initial_number_of_iterations = 200

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
            if tot_entropy_over_time[end] < 50 || tot_entropy_over_time[end-1] - tot_entropy_over_time[end] <= 0.5 || isnan(tot_entropy_over_time[end])
                break
            end
        end

        # number_of_iterations_of_ends = 100
        # rounds_for_ends = 2
        # variables_at_ends = [union(additional_variables, variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
        # factors_at_ends = [union(additional_factors, factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]
        # for i in 1:number_of_iterations_of_ends
        #     println(i)
        #     for var_name in variables_at_ends
        #         variable_to_factor_messages(variables[var_name], damping_factor)
        #     end
        #     for fact_name in factors_at_ends
        #         factor_to_variable_messages(factors[fact_name], damping_factor)
        #     end
        #     update_all_entropies(variables, variables_at_ends)
        #     push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))

        #     push!(tot_entropy_over_time, total_entropy_of_graph(variables))
        #     println(tot_entropy_over_time[end])
        #     if tot_entropy_over_time[end] < 100 || tot_entropy_over_time[end-1] - tot_entropy_over_time[end] <= 0.5
        #         break
        #     end
        # end

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
number_of_encryption_traces::Int64 = 1

num_successes = @distributed (+) for i in initial_key_number:final_key_number
    evaluate_key_number(number_of_bits, bits_per_template, dimensions_per_template, number_of_encryption_traces, i)
end

println(num_successes)