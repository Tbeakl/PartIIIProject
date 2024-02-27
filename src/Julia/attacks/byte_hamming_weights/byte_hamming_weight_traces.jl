include("../../belief_propagation/node.jl")
include("../../encryption/hamming_weight_probability_calculations.jl")

function byte_hamming_weight_value_to_function(hamming_position_table::Matrix{Bool})
    return function add_byte_hamming_weight_to_variable(value,
        variables::Dict{String,Variable{Factor}},
        factors::Dict{String,Factor{Variable}},
        bits_per_cluster::Int64,
        variable_and_count::String,
        run_number::Int64)
        clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
        for i in eachindex(value)
            hamming_value_likelihoods = likelihood_of_hamming_values_no_noise(8, bits_per_cluster, value[i])
            prob_dist_for_cluster = make_prob_distribution_from_hamming_likelihoods(hamming_value_likelihoods, hamming_position_table, bits_per_cluster)
            # prob_dist_for_cluster ./= sum(prob_dist_for_cluster)
            for j in 1:clusters_per_leakage_weight
                cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j, "_", run_number)
                cur_dist_name = string("f_", cur_var_name, "_dist")
                factors[cur_dist_name] = Factor{Variable}(cur_dist_name, LabelledArray(prob_dist_for_cluster, [cur_var_name]))
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