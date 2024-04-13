using HDF5, Base.Threads, StatsBase, Statistics, Distributions, LinearAlgebra
include("../attacks/byte_template_attacks/template_attack_traces.jl")


function get_prob_dist_of_vector(mean_vectors, noise, current_vector)
    likelihood_of_values = pdf(noise, mean_vectors' .- current_vector)
    return likelihood_of_values ./ sum(likelihood_of_values)
end

function noise_distribution_given_covaraince_matrix(scov)
    return MvNormal(Hermitian(scov))
end

function main()
    # Need to change these to work correctly for averaged versions of the templates
    fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_validation.hdf5", "r")

    all_intermediate_values = read(fid["intermediate_values"])
    downsampled_matrix = read(fid["downsampled_matrix"])

    number_of_bits_per_template = 8
    number_of_templates_per_intermediate_byte = 8 ÷ number_of_bits_per_template
    intermediate_value_index = 1933
    template_number = 1

    number_correct = 0
    println(template_number)
    template_path = string("D:\\ChaChaData\\attack_profiling\\initial_templates\\", intermediate_value_index, "_template.hdf5")
    intermediate_value_vector = (all_intermediate_values[:, intermediate_value_index] .>> (number_of_bits_per_template * (template_number - 1))) .& ((1 << number_of_bits_per_template) - 1)
    fid = h5open(template_path, "r")
    sample_bitmask = read(fid["downsampled_sample_bitmask"])
    template_projection = read(fid["projection"])
    mean_vectors = read(fid["class_means"])
    cov_matrix = read(fid["covariance_matrix"])
    close(fid)
    projected_vectors = downsampled_matrix[:, sample_bitmask] * template_projection
    noise_distribution = noise_distribution_given_covaraince_matrix(cov_matrix)
    values_predicted = zeros(Int64, 256)
    for j in eachindex(intermediate_value_vector)
        most_likely_value = findmax(get_prob_dist_of_vector(mean_vectors, noise_distribution, projected_vectors[j, :]))[2] - 1
        values_predicted[most_likely_value + 1] += 1
        if intermediate_value_vector[j] == most_likely_value
            number_correct += 1
        end
    end
    println(number_correct)
    # It is pretty random in the values it is predicting
    # p = scatter(size=(700, 700))
    # for i in 0:3
    #     scatter!(p, [projected_vectors[:, 1][intermediate_value_vector .== i]], [projected_vectors[:, 2][intermediate_value_vector .== i]], label=string("Samples ", i))
end
# scatter!(p, [mean_vectors[:, 1]], [mean_vectors[:, 2]])
main()