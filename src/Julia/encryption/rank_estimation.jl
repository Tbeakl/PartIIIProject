using DSP, StatsBase

include("../belief_propagation/node.jl")
include("../belief_propagation/messages.jl")

function make_log_likelihood_tables_for_key(variables::Dict{String,AbsVariable},
    bits_per_cluster::Int64)
    tables::Vector{Vector{Float64}} = []
    number_of_cluster_per_word = Int64(ceil(32 / bits_per_cluster))
    for i in 1:8
        for cluster_num in 1:number_of_cluster_per_word
            push!(tables, log2.(marginal(variables[string(i + 4, "_0_", cluster_num, "_1")])))
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

function chinese_remainder(inv_mods::Vector, Π::BigInt, divved_values::Vector, a::Vector)
    a = BigInt.(a)
    return mod(sum(a[i] * inv_mods[i] * divved_values[i] for i in eachindex(a)), Π)
end

function rank_estimate_key(key::Vector{UInt32},
    log_likelihood_tables::Vector{Vector{Float64}},
    number_of_bins::Int64,
    bits_per_cluster::Int64)
    min_likelihood = minimum(minimum.(log_likelihood_tables)) - 0.001
    max_likelihood = maximum(maximum.(log_likelihood_tables)) + 0.001
    log_likelihood_of_key = calculate_log_likelihood_of_key(key, log_likelihood_tables, bits_per_cluster)

    edges = range(min_likelihood, max_likelihood, number_of_bins + 1)
    # Some are missing from here like the top values if they are fractionally above
    histogram_bins = [(fit(Histogram, log_likelihood_tables[i], edges).weights) for i in 1:length(log_likelihood_tables)]
    # Need to work out how to do the convolution because it seems that the regular conv method does not work correctly
    # due to the large size of the integers therefore might need to use their suggested chinese remainder theorem method
    # and come up with a way for calculating the mid points of each of the bins, this would involve using the 20 largest prime numbers
    # below 10000 or potentially just look at the performance of using the ExactConvolution library for big integers 
    # because that might be fast enough BigInt

    # Need to perform for each of the primes the convolution mod that prime
    all_final_histograms::Vector{Vector{Int64}} = []
    largest_primes_below_10000 = [9811, 9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, 9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973]
    for prime in largest_primes_below_10000
        curr_hist_bins = histogram_bins[1]
        for i in 2:length(histogram_bins)
            curr_hist_bins = conv(histogram_bins[i], curr_hist_bins) .% prime
        end
        push!(all_final_histograms, curr_hist_bins)
    end
    histogram_matrix = permutedims(hcat(all_final_histograms...))
    largest_primes_below_10000 = BigInt.(largest_primes_below_10000)
    Π = prod(largest_primes_below_10000)
    divved_values = Π .÷ largest_primes_below_10000
    inv_mods = invmod.(divved_values, largest_primes_below_10000)

    # Need to use the chinese remainder theorem to calculate for each of the values of the final histogram the actual value they take
    actual_histogram_values = zeros(BigInt, size(histogram_matrix)[2])
    for i in eachindex(actual_histogram_values)
        actual_histogram_values[i] = chinese_remainder(inv_mods, Π, divved_values, histogram_matrix[:, i])
    end

    step_size = (max_likelihood - min_likelihood) / number_of_bins

    new_bin_mid_points = range(length(histogram_bins) * (min_likelihood + (step_size / 2)), length(histogram_bins) * (max_likelihood - (step_size / 2)), length(actual_histogram_values))
    # Not entirely sure what the correct solution is here for the number of bins which should be included because it is clearly not working
    # quite right I don't think but it potentially is
    # Because the upper end of the bins to the left are decreasing from the original upper end of the binds at the step rate as well 
    # so to get the overestimate can look to get the upper end of the bin included instead
    bin_upper_values = (-(length(actual_histogram_values) - 1):0 .* step_size) .+ (length(histogram_bins) * max_likelihood)
    # vector_of_bins_to_include = new_bin_mid_points .>= log_likelihood_of_key
    vector_of_bins_to_include = bin_upper_values .>= log_likelihood_of_key
    return sum(actual_histogram_values[vector_of_bins_to_include])
end

# Pretty sure this code is giving the wrong estiamtes out or there it is not using 
# enough bins or something to be useful but I think it is something else because there should be lots of
# keys in the top bins I think

# Pretty sure the bin midpoints are incorrect

likelihood_tables = make_log_likelihood_tables_for_key(variables, number_of_bits)
estimated_rank = rank_estimate_key(key, likelihood_tables, 500, number_of_bits)