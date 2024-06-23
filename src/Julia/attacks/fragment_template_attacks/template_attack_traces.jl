using Distributions, HDF5
include("../../belief_propagation/node.jl")
include("../../profiling/common_functions.jl")

function get_dist_of_vector(mean_vectors, noise, current_vector)
    likelihood_of_values = pdf(noise, mean_vectors' .- current_vector)
    if sum(likelihood_of_values) == 0
        return ones(size(mean_vectors)[1])
    end
    return likelihood_of_values
end

function get_log_likelihoods_dist_of_vector(mean_vectors, noise, current_vector)
    likelihood_of_values = logpdf(noise, mean_vectors' .- current_vector)
    return likelihood_of_values
end

function marginalise_prob_dist(original_dist::AbstractVector{Float64}, shift_amount_to_right::Int64, bits_in_marginalisation::Int64)
    output_dist = zeros(1 << bits_in_marginalisation)
    original_indexes_to_new = (((0:length(original_dist)-1) .>> shift_amount_to_right) .% (1 << bits_in_marginalisation)) .+ 1
    for i in 1:((1<<bits_in_marginalisation))
        output_dist[i] = sum(original_dist[original_indexes_to_new.==i])
    end
    return output_dist
end

function real_byte_template_path_to_function(base_path::String, bits_per_template::Int64, number_of_dimensions::Int64, traces)
    # Need mapping from the positions in the theoertical trace to template number
    val = 8 + 4 + 16 + 1 # Need to add on the base of the 
    trace_mapping = zeros(Int64, 1600)
    for i in 1:5:1600
        trace_mapping[i] = val
        val += 1
        trace_mapping[i+4] = val
        val += 1
    end
    templates_per_intermediate_value = 32 ÷ bits_per_template
    return function add_byte_template_to_variable(trace::Vector{Float32},
        position_in_trace::Int64,
        variables::Dict{String,AbsVariable},
        factors::Dict{String,AbsFactor},
        bits_per_cluster::Int64,
        variable_and_count::String,
        run_number::Int64)
        clusters_per_leakage_weight = Int64(ceil(bits_per_template / bits_per_cluster))
        for template_number in 1:templates_per_intermediate_value
            template_path = string(base_path, trace_mapping[position_in_trace], "_", template_number, "_template.hdf5")
            if ispath(template_path)
                fid = h5open(template_path, "r")
                projection_matrix = read(fid["projection"])[:, 1:number_of_dimensions]
                mean_vectors = read(fid["class_means"])[:, 1:number_of_dimensions]
                covaraince_matrix = read(fid["covariance_matrix"])[1:number_of_dimensions, 1:number_of_dimensions]
                # detailed_sample_bitmask = read(fid["detailed_sample_bitmask"])
                # sparse_sample_bitmask = read(fid["sparse_sample_bitmask"])
                sample_bitmask = read(fid["sample_bitmask"])
                close(fid)
                # noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
                noise = noise_distribution_given_covaraince_matrix(covaraince_matrix) #  ./ size(traces)[1]

                # position = mean(hcat(traces[:, detailed_sample_bitmask], traces[:, sparse_sample_bitmask]) * projection_matrix, dims=1)[1, :]
                position = mean(traces[:, sample_bitmask] * projection_matrix, dims=1)[1, :]
                prob_dist_for_template = get_dist_of_vector(mean_vectors, noise, position)
                # println(size(position))
                # prob_dist_for_template = zeros(1 << bits_per_template)
                # for cur_trace in eachrow(traces)
                #     values = vcat(cur_trace[detailed_sample_bitmask], cur_trace[sparse_sample_bitmask])
                #     prob_dist_for_template += get_log_likelihoods_dist_of_vector(mean_vectors, noise, (values'*projection_matrix)[1, :])
                # end
                # # First bring up the max value to zero then exponentiate it, probably need to have some dealing with the problems associated with 
                # # floating points and exp
                # prob_dist_for_template .-= maximum(prob_dist_for_template)
                # # Make all the values less than -700 equal to -700 to hopefully remove some floating point errors
                # prob_dist_for_template[prob_dist_for_template .< -50] .= -50
                # prob_dist_for_template = exp.(prob_dist_for_template)
                # prob_dist_for_template ./= sum(prob_dist_for_template)

                for j in 1:clusters_per_leakage_weight
                    cur_var_name = string(variable_and_count, "_", (template_number - 1) * clusters_per_leakage_weight + j, "_", run_number)
                    cur_dist_name = string("f_", cur_var_name, "_dist")
                    # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                    marginalised_dist = marginalise_prob_dist(prob_dist_for_template, (j - 1) * bits_per_cluster, bits_per_cluster)
                    factors[cur_dist_name] = Factor{AbsVariable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                    add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                    variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
                end
            end
        end
    end
end

function add_initial_key_distribution_from_leakage_traces_set_of_templates(
    variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    bits_per_cluster::Int64,
    base_path::String,
    bits_per_template::Int64,
    number_of_dimensions::Int64,
    traces::Matrix{Float32}
)
    templates_per_intermediate_value = 32 ÷ bits_per_template
    clusters_per_leakage_weight = Int64(ceil(bits_per_template / bits_per_cluster))
    for word_number in 5:12
        for template_number in 1:templates_per_intermediate_value
            fid = h5open(string(base_path, word_number - 4, "_", template_number, "_template.hdf5"), "r")
            covaraince_matrix = read(fid["covariance_matrix"])
            cur_dims = min(number_of_dimensions, size(covaraince_matrix)[1])
            covaraince_matrix = covaraince_matrix[1:cur_dims, 1:cur_dims]
            projection_matrix = read(fid["projection"])[:, 1:cur_dims]
            mean_vectors = read(fid["class_means"])[:, 1:cur_dims]
            # detailed_sample_bitmask = read(fid["detailed_sample_bitmask"])
            # sparse_sample_bitmask = read(fid["sparse_sample_bitmask"])
            sample_bitmask = read(fid["sample_bitmask"])
            close(fid)
            # noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
            noise = noise_distribution_given_covaraince_matrix(covaraince_matrix) #  ./ size(traces)[1]

            # position = mean(hcat(traces[:, detailed_sample_bitmask], traces[:, sparse_sample_bitmask]) * projection_matrix, dims=1)[1, :]
            position = mean(traces[:, sample_bitmask] * projection_matrix, dims=1)[1, :]
            prob_dist_for_template = get_dist_of_vector(mean_vectors, noise, position)

            # prob_dist_for_template = zeros(1 << bits_per_template)
            # for cur_trace in eachrow(traces)
            #     values = vcat(cur_trace[detailed_sample_bitmask], cur_trace[sparse_sample_bitmask])
            #     prob_dist_for_template += get_log_likelihoods_dist_of_vector(mean_vectors, noise, (values'*projection_matrix)[1, :])
            # end
            # # First bring up the max value to zero then exponentiate it, probably need to have some dealing with the problems associated with 
            # # floating points and exp
            # prob_dist_for_template .-= maximum(prob_dist_for_template)
            # prob_dist_for_template[prob_dist_for_template .< -50] .= -50
            # prob_dist_for_template = exp.(prob_dist_for_template)
            # prob_dist_for_template ./= sum(prob_dist_for_template)

            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(word_number, "_0_", (template_number - 1) * clusters_per_leakage_weight + j, "_1")
                cur_dist_name = string("f_", cur_var_name, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = marginalise_prob_dist(prob_dist_for_template, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{AbsVariable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function add_initial_nonce_distribution_from_leakage_traces_set_of_templates(
    variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    bits_per_cluster::Int64,
    base_path::String,
    bits_per_template::Int64,
    number_of_dimensions::Int64,
    traces::Matrix{Float32}
)
    templates_per_intermediate_value = 32 ÷ bits_per_template
    clusters_per_leakage_weight = Int64(ceil(bits_per_template / bits_per_cluster))
    for word_number in 14:16
        for template_number in 1:templates_per_intermediate_value
            fid = h5open(string(base_path, word_number - 4, "_", template_number, "_template.hdf5"), "r")
            covaraince_matrix = read(fid["covariance_matrix"])
            cur_dims = min(number_of_dimensions, size(covaraince_matrix)[1])
            covaraince_matrix = covaraince_matrix[1:cur_dims, 1:cur_dims]
            projection_matrix = read(fid["projection"])[:, 1:cur_dims]
            mean_vectors = read(fid["class_means"])[:, 1:cur_dims]
            # detailed_sample_bitmask = read(fid["detailed_sample_bitmask"])
            # sparse_sample_bitmask = read(fid["sparse_sample_bitmask"])
            sample_bitmask = read(fid["sample_bitmask"])
            close(fid)
            # noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
            noise = noise_distribution_given_covaraince_matrix(covaraince_matrix) #  ./ size(traces)[1]

            # position = mean(hcat(traces[:, detailed_sample_bitmask], traces[:, sparse_sample_bitmask]) * projection_matrix, dims=1)[1, :]
            position = mean(traces[:, sample_bitmask] * projection_matrix, dims=1)[1, :]
            prob_dist_for_template = get_dist_of_vector(mean_vectors, noise, position)

            # prob_dist_for_template = zeros(1 << bits_per_template)
            # for cur_trace in eachrow(traces)
            #     values = vcat(cur_trace[detailed_sample_bitmask], cur_trace[sparse_sample_bitmask])
            #     prob_dist_for_template += get_log_likelihoods_dist_of_vector(mean_vectors, noise, (values'*projection_matrix)[1, :])
            # end
            # # First bring up the max value to zero then exponentiate it, probably need to have some dealing with the problems associated with 
            # # floating points and exp
            # prob_dist_for_template .-= maximum(prob_dist_for_template)
            # prob_dist_for_template[prob_dist_for_template .< -50] .= -50
            # prob_dist_for_template = exp.(prob_dist_for_template)
            # prob_dist_for_template ./= sum(prob_dist_for_template)

            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(word_number, "_0_", (template_number - 1) * clusters_per_leakage_weight + j, "_1")
                cur_dist_name = string("f_", cur_var_name, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = marginalise_prob_dist(prob_dist_for_template, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{AbsVariable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function all_distribution_of_output(
    variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    bits_per_cluster::Integer,
    base_path::String,
    bits_per_template::Integer,
    number_of_dimensions::Integer,
    location_execution_counts::Vector{Int64},
    trace::Vector{Float32},
    encryption_run_number::Int64
)
    templates_per_intermediate_value = 32 ÷ bits_per_template
    clusters_per_leakage_weight = Int64(ceil(bits_per_template / bits_per_cluster))
    for word_number in 1:16
        for template_number in 1:templates_per_intermediate_value
            fid = h5open(string(base_path, 668 + word_number, "_", template_number, "_template.hdf5"), "r")
            covaraince_matrix = read(fid["covariance_matrix"])
            cur_dims = min(number_of_dimensions, size(covaraince_matrix)[1])
            covaraince_matrix = covaraince_matrix[1:cur_dims, 1:cur_dims]
            projection_matrix = read(fid["projection"])[:, 1:cur_dims]
            mean_vectors = read(fid["class_means"])[:, 1:cur_dims]
            sample_bitmask = read(fid["sample_bitmask"])
            close(fid)
            noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)

            position = (trace[sample_bitmask]'*projection_matrix)[1, :]
            prob_dist_for_template = get_dist_of_vector(mean_vectors, noise, position)

            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(word_number, "_", location_execution_counts[word_number], "_", (template_number - 1) * clusters_per_leakage_weight + j, "_", encryption_run_number)
                cur_dist_name = string("f_", cur_var_name, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = marginalise_prob_dist(prob_dist_for_template, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{AbsVariable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function add_initial_counter_distribution_from_leakage_traces_set_of_templates(
    variables::Dict{String,AbsVariable},
    factors::Dict{String,AbsFactor},
    bits_per_cluster::Int64,
    base_path::String,
    bits_per_template::Int64,
    number_of_dimensions::Int64,
    trace::Vector{Float32},
    encryption_run_number::Int64
)
    templates_per_intermediate_value = 32 ÷ bits_per_template
    clusters_per_leakage_weight = Int64(ceil(bits_per_template / bits_per_cluster))
    for template_number in 1:templates_per_intermediate_value
        fid = h5open(string(base_path, "9_", template_number, "_template.hdf5"), "r")
        covaraince_matrix = read(fid["covariance_matrix"])
        cur_dims = min(number_of_dimensions, size(covaraince_matrix)[1])
        covaraince_matrix = covaraince_matrix[1:cur_dims, 1:cur_dims]
        projection_matrix = read(fid["projection"])[:, 1:cur_dims]
        mean_vectors = read(fid["class_means"])[:, 1:cur_dims]
        sample_bitmask = read(fid["sample_bitmask"])
        close(fid)
        noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)

        position = (trace[sample_bitmask]'*projection_matrix)[1, :]
        prob_dist_for_template = get_dist_of_vector(mean_vectors, noise, position)

        for j in 1:clusters_per_leakage_weight
            cur_var_name = string("13_0_", (template_number - 1) * clusters_per_leakage_weight + j, "_", encryption_run_number)
            cur_dist_name = string("f_", cur_var_name, "_dist")
            # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
            marginalised_dist = marginalise_prob_dist(prob_dist_for_template, (j - 1) * bits_per_cluster, bits_per_cluster)
            factors[cur_dist_name] = Factor{AbsVariable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
            add_edge_between(variables[cur_var_name], factors[cur_dist_name])
            variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
        end
    end
end

function plot_distribution_of_values_and_means(
    base_path::String,
    number_of_dimensions::Int64,
    intermediate_value_index::Int64,
    template_number::Int64,
    value::Int64,
    traces,
    dim_1::Int64=1,
    dim_2::Int64=2
)
    fid = h5open(string(base_path, intermediate_value_index, "_", template_number, "_template.hdf5"), "r")
    covaraince_matrix = read(fid["covariance_matrix"])
    cur_dims = min(number_of_dimensions, size(covaraince_matrix)[1])
    covaraince_matrix = covaraince_matrix[1:cur_dims, 1:cur_dims]
    projection_matrix = read(fid["projection"])[:, 1:cur_dims]
    mean_vectors = read(fid["class_means"])[:, 1:cur_dims]
    # detailed_sample_bitmask = read(fid["detailed_sample_bitmask"])
    # sparse_sample_bitmask = read(fid["sparse_sample_bitmask"])
    sample_bitmask = read(fid["sample_bitmask"])
    close(fid)
    noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
    # current_value_matrix = hcat(traces[:, detailed_sample_bitmask], traces[:, sparse_sample_bitmask]) * projection_matrix
    current_value_matrix = traces[:, sample_bitmask] * projection_matrix
    position = mean(current_value_matrix, dims=1)[1, :]

    p = scatter(size=(1000, 1000))
    randomly_generated_values_around_mean = rand(noise, 10_000) .+ (mean_vectors[value+1, :])
    scatter!(p, [randomly_generated_values_around_mean[dim_1, :]], [randomly_generated_values_around_mean[dim_2, :]], label="Generated samples of particular value")
    scatter!(p, [mean_vectors[:, dim_1]], [mean_vectors[:, dim_2]], label="Class Means")
    scatter!(p, [current_value_matrix[:, dim_1]], [current_value_matrix[:, dim_2]], label="Actual samples of particular value")
    scatter!(p, [mean_vectors[value+1, dim_1]], [mean_vectors[value+1, dim_2]], markersize=10, label="Mean of class")
    scatter!(p, [position[dim_1]], [position[dim_2]], markersize=10, label="Mean of data")
    return p
end

function load_attack_trace(file_path::String, trace_number::Int64, encryption_run_number::Int64, path_to_mean::String)
    # Need to calculate the file in which it falls because each file contains 100
    # values with them and then also the encryption run number which needs to be done

    file_number = (trace_number ÷ 100)
    trace_number_in_file = (trace_number - 1) % 100
    clock_cycle_sample_number =  31 #46 #405
    number_of_samples_per_cycle = 100
    number_of_samples_to_average_over = 5

    fid = h5open(string(file_path, file_number, ".hdf5"), "r")
    base_trace_data = fid[string("power_", trace_number_in_file, "_", encryption_run_number)]
    key = UInt32.(read(base_trace_data["key"]))
    nonce = UInt32.(read(base_trace_data["nonce"]))
    counter = UInt32(read(base_trace_data["counter"])[1])

    raw_trace = read(base_trace_data)
    close(fid)
    # Need to downsample and align this trace to the mean
    fid = h5open(path_to_mean, "r")
    mean_trace = read(fid["mean_trace"])
    close(fid)

    trimmed_raw_trace = make_power_trace_trimmed_and_aligned_to_mean(mean_trace, raw_trace)
    trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(number_of_samples_per_cycle-clock_cycle_sample_number)-1)]
    downsampled_trace = Float32.(collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, number_of_samples_to_average_over))))
    return (key, nonce, counter, downsampled_trace)
end

function load_attack_trace_argmin_alignment(file_path::String, trace_number::Int64, encryption_run_number::Int64, path_to_mean::String)
    # Need to calculate the file in which it falls because each file contains 100
    # values with them and then also the encryption run number which needs to be done
    file_number = (trace_number ÷ 100)
    trace_number_in_file = (trace_number - 1) % 100
    clock_cycle_sample_number = 405
    number_of_samples_to_average_over = 50

    fid = h5open(string(file_path, file_number, ".hdf5"), "r")
    base_trace_data = fid[string("power_", trace_number_in_file, "_", encryption_run_number)]
    key = UInt32.(read(base_trace_data["key"]))
    nonce = UInt32.(read(base_trace_data["nonce"]))
    counter = UInt32(read(base_trace_data["counter"])[1])

    raw_trace = read(base_trace_data)
    close(fid)
    # Need to downsample and align this trace to the mean
    fid = h5open(path_to_mean, "r")
    mean_trace = read(fid["mean_trace"])
    mean_arg_min = argmin(mean_trace)
    close(fid)

    difference_between_mean_and_power = argmin(raw_trace) - mean_arg_min
    trimmed_raw_trace = raw_trace[50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
    trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(500-clock_cycle_sample_number)-1)]
    downsampled_trace = Float32.(collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, number_of_samples_to_average_over))))
    return (key, nonce, counter, downsampled_trace)
end