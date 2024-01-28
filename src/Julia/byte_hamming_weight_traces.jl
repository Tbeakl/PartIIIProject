using CryptoSideChannel
using Distributions
include("chacha.jl")
include("node.jl")
include("hamming_weight_probability_calculations.jl")

function mean_vector_for_value(value)
    output = zeros(Int64, 4)
    output[1] = Base.count_ones(value & 0xFF)
    output[2] = Base.count_ones((value >> 8) & 0xFF)
    output[3] = Base.count_ones((value >> 16) & 0xFF)
    output[4] = Base.count_ones((value >> 24) & 0xFF)
    return output
end

function encrypt_collect_trace(key::Vector{UInt32}, nonce::Vector{UInt32}, counter::UInt32)
    global trace
    trace = []
    closure = () -> trace

    key_logging = map(x -> Logging.SingleFunctionLog(x, closure, mean_vector_for_value), key)
    nonce_logging = map(x -> Logging.SingleFunctionLog(x, closure, mean_vector_for_value), nonce)
    counter_logging = Logging.SingleFunctionLog(counter, closure, mean_vector_for_value)

    encrypt(key_logging, nonce_logging, counter_logging)

    return copy(trace)
end

function add_dist_to_variable(values,
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64,
    variable_and_count::String,
    hamming_position_table::Matrix{Bool}
    )
    
    clusters_per_leakage_weight = Int64(ceil(8 / bits_per_cluster))
    for i in 1:length(values)
        hamming_value_likelihoods = likelihood_of_hamming_values_no_noise(8, bits_per_cluster, values[i])
        prob_dist_for_cluster = make_prob_distribution_from_hamming_likelihoods(hamming_value_likelihoods, hamming_position_table, bits_per_cluster)
        # prob_dist_for_cluster ./= sum(prob_dist_for_cluster)
        for j in 1:clusters_per_leakage_weight
            cur_var_name = string(variable_and_count, "_", (i - 1) * clusters_per_leakage_weight + j)
            cur_dist_name = string("f_", cur_var_name, "_dist")
            factors[cur_dist_name] = Factor(cur_dist_name, LabelledArray(prob_dist_for_cluster, [cur_var_name]))
            add_edge_between(variables[cur_var_name], factors[cur_dist_name])
            variables[cur_var_name].neighbour_index_to_avoid = length(variables[cur_var_name].neighbours)
        end
    end
end

# Currently need the clusters to fall exactly along with the leakages
function add_distribution_from_position_in_trace(trace::Vector{Any},
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    variable::Int64,
    hamming_position_table::Matrix{Bool}
    )
    global position_in_trace
    
    add_dist_to_variable(trace[position_in_trace], variables, factors, bits_per_cluster, string(variable, "_", location_execution_counts[variable]), hamming_position_table)
    # Need to pay attention to how to deal with rotations
    position_in_trace += 1
    location_execution_counts[variable] += 1
end

function add_qr_trace(trace::Vector{Any},
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64,
    a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64},
    hamming_position_table::Matrix{Bool}
    )
    global position_in_trace

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, a, hamming_position_table)
    # Need to do the part with shifting
    position_in_trace += 3
    if 16 % bits_per_cluster != 0
        # This means that the result of the rotation is directly going to be the output after the rotation
        location_execution_counts[d] += 1 
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, hamming_position_table)

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, c, hamming_position_table)
    position_in_trace += 3
    if 12 % bits_per_cluster != 0
        location_execution_counts[b] += 1
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, hamming_position_table)

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, a, hamming_position_table)
    # Need to do the part with shifting
    position_in_trace += 3
    if 8 % bits_per_cluster != 0
        # This means that the result of the rotation is directly going to be the output after the rotation
        location_execution_counts[d] += 1 
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, hamming_position_table)

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, c, hamming_position_table)
    position_in_trace += 3
    if 7 % bits_per_cluster != 0
        location_execution_counts[b] += 1
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, hamming_position_table)
end

function add_trace_to_factor_graph(trace::Vector{Any},
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64
    )
    global position_in_trace
    position_in_trace = 1
    location_execution_counts = ones(Int64, 16)
    hamming_position_table = table_for_hamming_values(bits_per_cluster)
    for i in 1:10
        add_qr_trace(trace, variables, factors, bits_per_cluster, 1, 5, 9, 13, location_execution_counts, hamming_position_table)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 2, 6, 10, 14, location_execution_counts, hamming_position_table)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 3, 7, 11, 15, location_execution_counts, hamming_position_table)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 4, 8, 12, 16, location_execution_counts, hamming_position_table)

        add_qr_trace(trace, variables, factors, bits_per_cluster, 1, 6, 11, 16, location_execution_counts, hamming_position_table)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 2, 7, 12, 13, location_execution_counts, hamming_position_table)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 3, 8, 9, 14, location_execution_counts, hamming_position_table)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 4, 5, 10, 15, location_execution_counts, hamming_position_table)
    end
    # Need to have another part for setting the values to be exactly what they were in the output of the function because
    # we assume that we have access to those actual values compared to just their leakage values
end

function add_initial_key_dist(variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64,
    key::Vector{UInt32})
    hamming_position_table = table_for_hamming_values(bits_per_cluster)
    for i in 1:8
        add_dist_to_variable(mean_vector_for_value(key[i]), variables, factors, bits_per_cluster, string(i + 4, "_", 0), hamming_position_table)
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