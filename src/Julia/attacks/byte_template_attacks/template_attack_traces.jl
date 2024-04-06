using Distributions, HDF5
include("../../belief_propagation/node.jl")

function generate_mean_vectors(noise::Distribution, signal_to_noise_ratio::Real, max_value::Int64)
    return signal_to_noise_ratio .* rand(noise, max_value + 1)
end

function generate_mean_vectors_based_on_hamming_weights(distribution_from_weight::Distribution, number_of_bits::Int64)
    # The basic plan with this is to have a mean vector associated with each bit being set and then use the digits to add
    # these together to make the basic mean vector for different values
    bit_mean_vectors = rand(distribution_from_weight, number_of_bits) .+ 1
    binary_representation = digits.((0:(1<<number_of_bits)-1), base=2, pad=number_of_bits)
    output = zeros(size(bit_mean_vectors)[1], 1 << number_of_bits)
    for i in eachindex(binary_representation)
        output[:, i] = sum(bit_mean_vectors[:, Bool.(binary_representation[i])], dims=2)
    end
    return output
end

function put_value_into_noisy_space(mean_vectors::AbstractMatrix{Float64}, noise::Distribution, value::Int64)
    return mean_vectors[:, value+1] .+ rand(noise)
end

function get_prob_dist_of_vector(mean_vectors, noise, current_vector)
    likelihood_of_values = pdf(noise, mean_vectors' .- current_vector)
    return likelihood_of_values ./ sum(likelihood_of_values)
end

function make_prob_dist_for_byte(mean_vectors::AbstractMatrix{Float64}, noise::Distribution, value::Int64)
    likelihood_of_values = pdf(noise, mean_vectors .- put_value_into_noisy_space(mean_vectors, noise, value))
    return likelihood_of_values ./ sum(likelihood_of_values)
end

function marginalise_prob_dist(original_dist::AbstractVector{Float64}, shift_amount_to_right::Int64, bits_in_marginalisation::Int64)
    output_dist = zeros(1 << bits_in_marginalisation)
    original_indexes_to_new = (((0:length(original_dist)-1) .>> shift_amount_to_right) .% (1 << bits_in_marginalisation)) .+ 1
    for i in 1:((1<<bits_in_marginalisation))
        output_dist[i] = sum(original_dist[original_indexes_to_new.==i])
    end
    return output_dist
end

function byte_template_value_to_function(mean_vectors::AbstractMatrix{Float64}, noise::Distribution)
    return function add_byte_template_to_variable(value,
        variables::Dict{String,Variable{Factor}},
        factors::Dict{String,Factor{Variable}},
        bits_per_cluster::Int64,
        variable_and_count::String,
        run_number::Int64,
        version_to_add_include::Int64=run_number)

        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in 1:4
            prob_dist_for_byte = make_prob_dist_for_byte(mean_vectors, noise, value[i])
            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j, "_", version_to_add_include)
                cur_dist_name = string("f_", variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j, "_", run_number, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = marginalise_prob_dist(prob_dist_for_byte, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function byte_template_path_to_function(base_path::String, noise::Distribution, number_of_templates::Int64)
    return function add_byte_template_to_variable(value,
        variables::Dict{String,Variable{Factor}},
        factors::Dict{String,Factor{Variable}},
        bits_per_cluster::Int64,
        variable_and_count::String,
        run_number::Int64,
        version_to_add_include::Int64=run_number)
        mean_vectors = transpose(npzread(string(base_path, lpad(string(Base.rand(0:(number_of_templates-1))), 3, "0"), ".npy")))
        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in 1:4
            prob_dist_for_byte = make_prob_dist_for_byte(mean_vectors, noise, value[i])
            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j, "_", version_to_add_include)
                cur_dist_name = string("f_", variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j, "_", run_number, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = marginalise_prob_dist(prob_dist_for_byte, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function real_byte_template_path_to_function(base_path::String)
    # Need mapping from the positions in the theoertical trace to template number
    val = 1 + 112 # Need to add on the base of the 
    trace_mapping = zeros(Int64, 1600)
    for i in 1:5:1600
        trace_mapping[i] = val
        val += 4
        trace_mapping[i+4] = val
        val += 4
    end
    return function add_byte_template_to_variable(trace::Vector{Float32},
        position_in_trace::Int64,
        variables::Dict{String,Variable{Factor}},
        factors::Dict{String,Factor{Variable}},
        bits_per_cluster::Int64,
        variable_and_count::String,
        run_number::Int64)
        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in 1:4
            fid = h5open(string(base_path, trace_mapping[position_in_trace] + i - 1, "_template.hdf5"), "r")
            projection_matrix = read(fid["projection"])
            mean_vectors = read(fid["class_means"])
            covaraince_matrix = read(fid["covariance_matrix"])
            sample_bitmask = read(fid["downsampled_sample_bitmask"])
            close(fid)
            noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)

            # println("Mean Vectors")
            # println(size(mean_vectors))
            # println("Projected Values")
            # println(size(trace[sample_bitmask]' * projection_matrix))

            prob_dist_for_byte = get_prob_dist_of_vector(mean_vectors, noise, (trace[sample_bitmask]'*projection_matrix)[1, :])
            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j, "_", run_number)
                cur_dist_name = string("f_", cur_var_name, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = marginalise_prob_dist(prob_dist_for_byte, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function add_initial_key_distribution_from_leakage_trace(trace::Vector{Float32},
    variables::Dict{String,Variable{Factor}},
    factors::Dict{String,Factor{Variable}},
    bits_per_cluster::Int64,
    run_number::Int64,
    base_path::String
    )
    clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
    for word_number in 5:12
        for i in 1:4
            for j in 1:clusters_per_leakage_weight
                fid = h5open(string(base_path, (word_number - 5) * 4 + i, "_", j, "_template.hdf5"), "r")
                projection_matrix = read(fid["projection"])
                mean_vectors = read(fid["class_means"])
                covaraince_matrix = read(fid["covariance_matrix"])
                sample_bitmask = read(fid["downsampled_sample_bitmask"])
                close(fid)
                noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
                prob_dist_for_byte = get_prob_dist_of_vector(mean_vectors, noise, (trace[sample_bitmask]'*projection_matrix)[1, :])

                cur_var_name = string(word_number, "_0_", (i - 1) * clusters_per_leakage_weight + j, "_", run_number)
                cur_dist_name = string("f_", cur_var_name, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = prob_dist_for_byte # marginalise_prob_dist(prob_dist_for_byte, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function add_initial_key_distribution_from_simulated_leakage(key_values::Vector,
    variables::Dict{String,Variable{Factor}},
    factors::Dict{String,Factor{Variable}},
    bits_per_cluster::Int64,
    run_number::Int64,
    base_path::String,
    number_of_reads::Int64
    )
    clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
    for word_number in 5:12
        for i in 1:4
            for j in 1:clusters_per_leakage_weight
                fid = h5open(string(base_path, (word_number - 5) * 4 + i, "_", j, "_template.hdf5"), "r")
                mean_vectors = read(fid["class_means"])
                covaraince_matrix = read(fid["covariance_matrix"])
                close(fid)
                noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
                value = (key_values[word_number - 4][i] .>> (bits_per_cluster * (j - 1))) .& ((1 << bits_per_cluster) - 1)
                leakage_output = mean(rand(noise, number_of_reads) .+ mean_vectors[value + 1, :], dims=2)
                noise_of_means = noise_distribution_given_covaraince_matrix(covaraince_matrix ./ number_of_reads)
                prob_dist_for_byte = get_prob_dist_of_vector(mean_vectors, noise_of_means, leakage_output)
                cur_var_name = string(word_number, "_0_", (i - 1) * clusters_per_leakage_weight + j, "_", run_number)
                cur_dist_name = string("f_", cur_var_name, "_dist")
                # Marginalise out the prob dist for this particular cluster, where cluster 1 is the LSB
                marginalised_dist = prob_dist_for_byte #marginalise_prob_dist(prob_dist_for_byte, (j - 1) * bits_per_cluster, bits_per_cluster)
                factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(marginalised_dist, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

function load_attack_trace(trace_number::Int64, encryption_run_number::Int64)
    # Need to calculate the file in which it falls because each file contains 100
    # values with them and then also the encryption run number which needs to be done

    file_number = (trace_number ÷ 100) + 1
    trace_number_in_file = (trace_number - 1) % 100
    clock_cycle_sample_number = 405
    number_of_samples_to_average_over = 10

    fid = h5open(string("D:\\ChaChaData\\captures\\ChaChaRecordings\\recording_attack_", file_number, ".hdf5"), "r")
    base_trace_data = fid[string("power_", trace_number_in_file, "_", encryption_run_number)]
    key = UInt32.(read(base_trace_data["key"]))
    nonce = UInt32.(read(base_trace_data["nonce"]))
    counter = UInt32(read(base_trace_data["counter"])[1])

    raw_trace = read(base_trace_data)
    close(fid)
    # Need to downsample and align this trace to the mean
    fid = h5open("D:\\ChaChaData\\attack_profiling\\mean_trace.hdf5", "r")
    mean_trace = read(fid["mean_trace"])
    mean_arg_min = argmin(mean_trace)
    close(fid)

    difference_between_mean_and_power = argmin(raw_trace) - mean_arg_min
    trimmed_raw_trace = raw_trace[50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
    trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(500-clock_cycle_sample_number)-1)]
    downsampled_trace = Float32.(collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, number_of_samples_to_average_over))))
    return (key, nonce, counter, downsampled_trace)
end