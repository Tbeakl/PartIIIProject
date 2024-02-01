using DSP, StatsBase, ExactConvolution

include("../belief_propagation/node.jl")
include("../belief_propagation/messages.jl")

function make_log_likelihood_tables_for_key(variables::Dict{String,Variable{Factor}},
    bits_per_cluster::Int64)
    tables::Vector{Vector{Float64}} = []
    number_of_cluster_per_word = Int64(ceil(32 / bits_per_cluster))
    for i in 1:8
        for cluster_num in 1:number_of_cluster_per_word
            push!(tables, log2.(marginal(variables[string(i + 4, "_0_", cluster_num)])))
        end
    end
    return tables
end

function calculate_log_likelihood_of_key(key::Vector{UInt32},
    log_likelihood_tables::Vector{Vector{Float64}},
    bits_per_cluster::Int64)
    total_log_likelihood = 0

    clusters_per_key_word = Int64(ceil(32 / bits_per_cluster))
    for i in 1:8
        for cluster_num in 1:clusters_per_key_word
            total_log_likelihood += log_likelihood_tables[clusters_per_key_word*(i-1)+cluster_num][extract_cluster_value_from_word(key[i], bits_per_cluster, cluster_num)+1]
        end
    end
    return total_log_likelihood
end

function extract_cluster_value_from_word(value::UInt32,
    bits_per_cluster::Int64,
    cluster_number::Int64)
    return (value >> (bits_per_cluster * (cluster_number - 1))) & ((1 << bits_per_cluster) - 1)
end

function rank_estimate_key(key::Vector{UInt32},
    log_likelihood_tables::Vector{Vector{Float64}},
    number_of_bins::Int64,
    bits_per_cluster::Int64)
    min_likelihood = minimum(minimum.(log_likelihood_tables))
    max_likelihood = maximum(maximum.(log_likelihood_tables))
    log_likelihood_of_key = calculate_log_likelihood_of_key(key, log_likelihood_tables, bits_per_cluster)

    edges = range(min_likelihood, max_likelihood, number_of_bins + 1)
    histogram_bins = [BigInt.(fit(Histogram, log_likelihood_tables[i], edges).weights) for i in 1:length(log_likelihood_tables)]
    # Need to work out how to do the convolution because it seems that the regular conv method does not work correctly
    # due to the large size of the integers therefore might need to use their suggested chinese remainder theorem method
    # and come up with a way for calculating the mid points of each of the bins, this would involve using the 20 largest prime numbers
    # below 10000 or potentially just look at the performance of using the ExactConvolution library for big integers 
    # because that might be fast enough
    # largest_primes_below_10000 = [9811,9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973]

    # Need to calculate the mid point likelihoods (might be that the bins are reversed in order compared to what I think they should be
    # because currently the top bins have a lot of weight which is completely infeasible due to no actual key going into those bins)
    curr_hist_bins = histogram_bins[1]
    for i in 2:length(histogram_bins)
        curr_hist_bins = exact_conv(BigInt, histogram_bins[i], curr_hist_bins)
    end

    step_size = (max_likelihood - min_likelihood) / number_of_bins
    
    # Not entirely sure if this is the correct method for calculating the mid points of the histogram bins
    new_bin_mid_points = range(length(histogram_bins) * (min_likelihood + (step_size / 2)), length(histogram_bins) * (max_likelihood - (step_size / 2)), length(curr_hist_bins))
    vector_of_bins_to_include = new_bin_mid_points .>= (log_likelihood_of_key - (step_size / 2)) # Also may want to change this
    return sum(curr_hist_bins[vector_of_bins_to_include])
end


likelihood_tables = make_log_likelihood_tables_for_key(variables, number_of_bits)
estimated_rank = rank_estimate_key(key, likelihood_tables, 200, 2)