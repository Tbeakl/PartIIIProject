using CryptoSideChannel
using Distributions
include("chacha.jl")
include("node.jl")
include("hamming_weight_probability_calculations.jl")

function noisy_byte_hamming_weight_value_to_function(hamming_position_table::Matrix{Bool}, noise_distribution::Distribution)
    return function add_noisy_byte_hamming_weight_to_variable(value,
        variables::Dict{String,Variable},
        factors::Dict{String,Factor},
        bits_per_cluster::Int64,
        variable_and_count::String)

        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in 1:length(value)
            hamming_value_likelihoods = likelihoods_of_hamming_values(Normal(0., sqrt(noise_distribution.Σ[i,i])), 8, bits_per_cluster, value[i])
            prob_dist_for_cluster = make_prob_distribution_from_hamming_likelihoods(hamming_value_likelihoods, hamming_position_table, bits_per_cluster)
            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j)
                cur_dist_name = string("f_", cur_var_name, "_dist")
                factors[cur_dist_name] = Factor(cur_dist_name, LabelledArray(prob_dist_for_cluster, [cur_var_name]))
                add_edge_between(variables[cur_var_name], factors[cur_dist_name])
                variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
            end
        end
    end
end

# How the full trace (1616) is broken down

# The first 1600 are from the 20 rounds, so that is 80 per round and 20 per quarter round

# Inside a quarter round is broken into 4 parts
# The add has a single part in the trace
# The exclusive-or also has a single part in the trace
# The rotation has three in the trace but I think that only the final should be used
# because that is more fundamental to the actual algorithm

# Last 16 elements from the addition afterwards