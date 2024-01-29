using Distributions
include("../../belief_propagation/node.jl")

function generate_mean_vectors(noise::Distribution, signal_to_noise_ratio::Real, max_value::Int64)
    return signal_to_noise_ratio .* rand(noise, max_value + 1)
end

function put_value_into_noisy_space(mean_vectors::Matrix{Float64}, noise::Distribution, value::Int64)
    return mean_vectors[:, value + 1] .+ rand(noise)
end

function make_prob_dist_for_byte(mean_vectors::Matrix{Float64}, noise::Distribution, value::Int64)
    likelihood_of_values = pdf(noise, mean_vectors .- put_value_into_noisy_space(mean_vectors, noise, value))
    return likelihood_of_values ./ sum(likelihood_of_values)
end

function marginalise_prob_dist(original_dist::Vector{Float64}, shift_amount_to_right::Int64, bits_in_marginalisation::Int64)
    output_dist = zeros(1 << bits_in_marginalisation)
    original_indexes_to_new = (((0:length(original_dist) - 1) .>> shift_amount_to_right) .% (1 << bits_in_marginalisation)) .+ 1
    for i in 1:((1 << bits_in_marginalisation))
        output_dist[i] = sum(original_dist[original_indexes_to_new .== i])
    end
    return output_dist
end

function byte_template_value_to_function(mean_vectors::Matrix{Float64}, noise::Distribution)
    return function add_byte_template_to_variable(value,
        variables::Dict{String, Variable{Factor}},
        factors::Dict{String, Factor{Variable}},
        bits_per_cluster::Int64,
        variable_and_count::String)
        
        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in 1:4
            prob_dist_for_byte = make_prob_dist_for_byte(mean_vectors, noise, value[i])
            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j)
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

function add_key_dist_byte_templates(variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    bits_per_cluster::Int64,
    key::Vector{UInt32},
    add_byte_template_to_variable)
    key_leakage = byte_values_for_input.(key)
    for i in 1:8
        add_byte_template_to_variable(key_leakage[i], variables, factors, bits_per_cluster, string(i + 4, "_", 0))
    end
end
