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
    signal_to_noise_ratio::Float64)

    path_to_data = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/"
    base_evaluation_of_ranks = string(path_to_data, "evaluation/simulation_unknown_output_nonce_counter/", signal_to_noise_ratio, "/")
    mkpath(base_evaluation_of_ranks)
    damping_factor::Float64 = 0.95

    # Load all the keys etc.
    all_keys::Vector{Vector{UInt32}} = []
    all_nonces::Vector{Vector{UInt32}} = []
    all_counters::Vector{UInt32} = []

    keyset_fid = h5open(string(path_to_data, "evaluation/key_sets/100_random_keys_nonces_counter.hdf5"), "r")
    for i in 1:100
        push!(all_keys, read(keyset_fid[string("key_", i)]))
        push!(all_nonces, read(keyset_fid[string("nonce_", i)]))
        push!(all_counters, read(keyset_fid[string("counter_", i)]))
    end
    close(keyset_fid)

    for key_number in 1:100
        Random.seed!(1234)

        noise = noise_distribution_fixed_standard_dev(1.0, dimensions)
        mean_vectors = generate_mean_vectors(noise, signal_to_noise_ratio, 255)
        add_byte_key_template_to_variable = byte_template_value_to_function(mean_vectors, noise)
        add_byte_intermediate_add_template_to_variable = byte_template_value_to_function(mean_vectors, noise)
        add_byte_intermediate_rot_template_to_variable = byte_template_value_to_function(mean_vectors, noise)

        result_path = string(base_evaluation_of_ranks, key_number, ".hdf5")
        if !ispath(result_path)
            variables = Dict{String,AbsVariable}()
            factors = Dict{String,AbsFactor}()
            variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
            factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
            adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
            println("Key number: ", key_number)
            encryption_trace = encrypt_collect_trace(key, nonce, counter, byte_values_for_input)
            encryption_output = encrypt(key, nonce, counter)
            location_execution_counts = zeros(Int64, 16)

            chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)
            add_starting_constant_values(variables, factors, number_of_bits, 1)

            println("Adding nonce, counter and key dist")
            add_initial_nonce_and_counter_dist(variables, factors, number_of_bits, byte_values_for_input.(nonce), byte_values_for_input(counter), 1, add_byte_key_template_to_variable)
            add_initial_key_dist(variables, factors, number_of_bits, byte_values_for_input.(key), 1, add_byte_key_template_to_variable)
            println("Added nonce, counter and key dist")

            println("Starting to add the factors for the trace")
            add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, 1, add_byte_intermediate_add_template_to_variable, add_byte_intermediate_rot_template_to_variable)
            println("Added the factors for the trace")

            # Also want to add in the known outputs of the encryption as part of the model
            println("Starting adding output distribution")
            for i in 1:16
                add_byte_key_template_to_variable(byte_values_for_input(encryption_output[i]), variables, factors, number_of_bits, string(i, "_", location_execution_counts[i]), 1)
            end
            println("Added output distribution")

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
                if tot_entropy_over_time[end] < 10 || tot_entropy_over_time[end-1] - tot_entropy_over_time[end] <= 0.5
                    break
                end
            end

            final_likelihood_tables = make_log_likelihood_tables_for_key(variables, number_of_bits)
            # Replace all the -Inf with like -800 because of the issues associated with expoentiatation
            final_estimated_rank = rank_estimate_key(key, final_likelihood_tables, 500, number_of_bits)

            # read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1) for i in 1:8]
            # for i in 1:8
            #     println(string(read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1), base=16))
            # end

            result_fid = h5open(result_path, "w")
            result_fid["initial_likelihood_tables"] = reduce(vcat, transpose.(initial_likelihood_tables))
            initial_estimated_rank = max(1, initial_estimated_rank)
            result_fid["initial_estimated_rank_log2"] = Float64(log2(initial_estimated_rank))
            result_fid["final_likelihood_tables"] = reduce(vcat, transpose.(final_likelihood_tables))
            final_estimated_rank = max(1, final_estimated_rank)
            result_fid["final_estimated_rank_log2"] = Float64(log2(final_estimated_rank))
            result_fid["entropy_over_time"] = tot_entropy_over_time
            close(result_fid)
        end
    end
    1
end

number_of_bits::Int64 = 8
num_successes = @distributed (+) for i in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
    evaluate_key_number(number_of_bits, i)
end

println(num_successes)