using Distributions
include("../../belief_propagation/node.jl")

function generate_mean_vectors(noise::Distribution, signal_to_noise_ratio::Real, max_value::Int64)
    return signal_to_noise_ratio .* rand(noise, max_value + 1)
end

function generate_mean_vectors_based_on_hamming_weights(distribution_from_weight::Distribution, number_of_bits::Int64)
    # The basic plan with this is to have a mean vector associated with each bit being set and then use the digits to add
    # these together to make the basic mean vector for different values
    bit_mean_vectors = rand(distribution_from_weight, number_of_bits) .+ 1
    binary_representation = digits.((0:(1 << number_of_bits) - 1), base = 2, pad = number_of_bits)
    output = zeros(size(bit_mean_vectors)[1], 1 << number_of_bits)
    for i in eachindex(binary_representation)
        output[:, i] = sum(bit_mean_vectors[:, Bool.(binary_representation[i])], dims=2)
    end
    return output
end

function put_value_into_noisy_space(mean_vectors::AbstractMatrix{Float64}, noise::Distribution, value::Int64)
    return mean_vectors[:, value + 1] .+ rand(noise)
end

function make_prob_dist_for_byte(mean_vectors::AbstractMatrix{Float64}, noise::Distribution, value::Int64)
    likelihood_of_values = pdf(noise, mean_vectors .- put_value_into_noisy_space(mean_vectors, noise, value))
    return likelihood_of_values ./ sum(likelihood_of_values)
end

function marginalise_prob_dist(original_dist::AbstractVector{Float64}, shift_amount_to_right::Int64, bits_in_marginalisation::Int64)
    output_dist = zeros(1 << bits_in_marginalisation)
    original_indexes_to_new = (((0:length(original_dist) - 1) .>> shift_amount_to_right) .% (1 << bits_in_marginalisation)) .+ 1
    for i in 1:((1 << bits_in_marginalisation))
        output_dist[i] = sum(original_dist[original_indexes_to_new .== i])
    end
    return output_dist
end

function byte_template_value_to_function(mean_vectors::AbstractMatrix{Float64}, noise::Distribution)
    return function add_byte_template_to_variable(value,
        variables::Dict{String, Variable{Factor}},
        factors::Dict{String, Factor{Variable}},
        bits_per_cluster::Int64,
        variable_and_count::String,
        run_number::Int64)

        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in 1:4
            prob_dist_for_byte = make_prob_dist_for_byte(mean_vectors, noise, value[i])
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
