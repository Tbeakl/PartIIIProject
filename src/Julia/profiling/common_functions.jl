using Statistics, StatsBase, LinearAlgebra, MultivariateStats

function dilate_infront(bit_vector::AbstractVector, dilation_amount::Int64)
    if dilation_amount == 0
        return bit_vector
    end
    bit_vector = copy(bit_vector)
    end_locations_of_dilations = findall(x -> x == 1, diff(bit_vector))
    for loc in end_locations_of_dilations
        bit_vector[max(loc - dilation_amount + 1, 1):loc] .= true
    end
    return bit_vector
end

function dilate_after(bit_vector::AbstractVector, dilation_amount::Int64)
    if dilation_amount == 0
        return bit_vector
    end
    bit_vector = copy(bit_vector)
    start_locations_of_dilations = findall(x -> x == -1, diff(bit_vector))
    end_index = length(bit_vector)
    for loc in start_locations_of_dilations
        bit_vector[loc:min(loc + dilation_amount, end_index)] .= true
    end
    return bit_vector
end


function calculate_within_class_scatter(data_samples, labels)
    # Data samples has the dimensions of (n, d) where n is the number of samples and d is the dimensions of the data
    output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
    for i in unique(labels)
        output_matrix += StatsBase.cov(data_samples[labels.==i, :], corrected=true)
    end
    return output_matrix ./ length(unique(labels))
end

# function calculate_within_class_scatter(data_samples, labels, means)
#     # Data samples has the dimensions of (n, d) where n is the number of samples and d is the dimensions of the data
#     output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
#     for i in unique(labels)
#         part_of_interest = data_samples[labels.==i, :] .- means[i+1, :]'
#         output_matrix += Statistics.cov(StatsBase.SimpleCovariance(; corrected=true), part_of_interest; mean=zeros(size(part_of_interest)))
#     end
#     return output_matrix ./ length(unique(labels))
# end

# function calculate_between_class_scatter(data_samples, labels, means)
#     # Data samples has the dimensions of (d, n) where n is the number of samples and d is the dimensions of the data
#     output_matrix = zeros(size(data_samples)[1], size(data_samples)[1])
#     overall_mean = mean(means, dims=2)[:, 1]
#     for i in unique(labels)
#         vals_to_include = labels .== i
#         cur_mean = means[:, i+1]
#         output_matrix += ((cur_mean - overall_mean) * (cur_mean - overall_mean)') .* sum(vals_to_include)
#     end
#     return output_matrix ./ length(unique(labels))
# end

function calculate_within_class_scatter(data_samples, labels, means)
    # Data samples has the dimensions of (n, d) where n is the number of samples and d is the dimensions of the data
    output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
    @assert size(data_samples)[1] == length(labels)
    @assert maximum(labels) < size(means)[1]
    @assert minimum(labels) >= 0
    @inbounds for i in eachindex(labels)
        @inbounds part_of_interest = data_samples[i, :] .- means[labels[i]+1, :]'
        output_matrix += part_of_interest * part_of_interest' #Statistics.cov(StatsBase.SimpleCovariance(; corrected=true), part_of_interest; mean=zeros(size(part_of_interest)))
    end
    return output_matrix ./ length(labels)
end

function alternative_within_class_scatter(data_samples, labels, means)
    @assert size(data_samples)[2] == length(labels)
    @assert maximum(labels) < size(means)[1]
    @assert minimum(labels) >= 0
    @inbounds for i in eachindex(labels)
        data_samples[:, i] -= means[labels[i]+1, :]
    end
    return (data_samples * data_samples') ./ length(labels)
end

function calculate_between_class_scatter(labels, means)
    output_matrix = zeros(size(means)[2], size(means)[2])
    overall_mean = mean(original_mean_vectors, dims=1)[1, :]
    for i in 0:size(means)[1]
        vals_to_include = labels .== i
        cur_mean = means[i+1, :]
        output_matrix += sum(vals_to_include) .* ((cur_mean - overall_mean) * (cur_mean - overall_mean)')
    end
    return output_matrix ./ length(labels)
end


function alternative_calculate_between_class_scatter(labels, means)
    overall_mean = mean(means, dims=1)[1, :]
    # label_counts = zeros(Int64, size(means)[1])
    matrix_of_vectors = zeros(size(means)[2], size(means)[1])
    for i in 0:(size(means)[1]-1)
        # label_counts[i + 1] = sum(labels .== i)
        matrix_of_vectors[:, i+1] = sum(labels .== i) .* (means[i+1, :] - overall_mean)
    end
    # for i in 0:size(means)[1]
    #     vals_to_include = labels .== i
    #     cur_mean = means[i + 1, :]
    #     output_matrix += sum(vals_to_include) .* ((cur_mean - overall_mean) * (cur_mean - overall_mean)')
    # end
    return (matrix_of_vectors * matrix_of_vectors') ./ length(labels)
end

# regularize a symmetric matrix
function regularize_symmat!(A::AbstractMatrix{T}, lambda::Real) where {T<:Real}
    if lambda > 0
        emax = eigmax(Symmetric(A))
        add_diag!(A, emax * lambda)
    end
    return A
end

function add_diag!(A::AbstractMatrix, v::Real)
    # add v to diagonal of A
    m = size(A, 1)
    n = size(A, 2)
    @assert m == n
    if v != zero(v)
        for i = 1:n
            @inbounds A[i, i] += v
        end
    end
    return A
end

function make_intermediate_value_matrix(datasets::AbstractVector, index::Integer)
    return reduce(vcat, [dset[:, index] for dset in datasets])
end

function make_dataset_of_values_matrix(datasets::AbstractVector, bitmask::BitVector)
    return reduce(vcat, [dset[:, bitmask] for dset in datasets])
end

function make_power_trace_trimmed_and_aligned_to_mean(mean_trace::Vector{Float32}, power_trace::Vector{Int16})
    base_difference_between_mean_and_power = argmin(power_trace[begin+50:end-50]) - argmin(mean_trace)
    lags_to_try = (-5:5) .+ base_difference_between_mean_and_power
    difference_between_mean_and_power = lags_to_try[argmax(crosscor(mean_trace, power_trace[begin+50:end-50], lags_to_try))]
    return power_trace[begin+50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
end

function get_prob_dist_of_vector(mean_vectors, noise, current_vector)
    likelihood_of_values = logpdf(noise, mean_vectors' .- current_vector)
    return likelihood_of_values# ./ sum(likelihood_of_values)
end

function get_number_of_guesses(mean_vectors, noise, value, projected_vector)
    permutation_of_results = sortperm(logpdf(noise, mean_vectors .- projected_vector); rev=true)
    return findfirst(x -> x == value + 1, permutation_of_results)
end

function noise_distribution_given_covaraince_matrix(scov)
    return MvNormal(Hermitian(scov))
end

function noise_distribution_fixed_standard_dev(standard_deviation, dimensions)
    variance = standard_deviation * standard_deviation
    covariance_matrix = ones(dimensions) .* variance
    return MvNormal(covariance_matrix)
end

function get_success_rate(means, noise, values, projected_vectors)
    number_of_successes = 0
    for i in eachindex(values)
        most_likely_value = findmax(get_prob_dist_of_vector(means, noise, projected_vectors[i, :]))[2] - 1
        if values[i] == most_likely_value
            number_of_successes += 1
        end
    end
    return number_of_successes / length(values)
end

function get_guessing_entropy(means, noise, values, projected_vectors)
    guessing_entropy = sum(get_number_of_guesses.(Ref(means'), Ref(noise), values, eachrow(projected_vectors)))
    return guessing_entropy / length(values)
end
