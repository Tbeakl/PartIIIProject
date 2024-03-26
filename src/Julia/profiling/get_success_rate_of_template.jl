using HDF5, Base.Threads, StatsBase, Statistics
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")


fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_validation.hdf5", "r")

all_intermediate_values = read(fid["intermediate_values"])
downsampled_matrix = read(fid["downsampled_matrix"])
number_of_templates = 2684

all_correct_counts = zeros(Int64, number_of_templates)
template_number = 1000
# for template_number in 1:number_of_templates
    println(template_number)
    template_path = string("D:\\ChaChaData\\attack_profiling\\initial_templates\\", template_number, "_template.hdf5")
    intermediate_value_vector = all_intermediate_values[:, template_number]
    fid = h5open(template_path, "r")
    sample_bitmask = read(fid["downsampled_sample_bitmask"])
    template_projection = read(fid["projection"])
    mean_vectors = read(fid["class_means"])
    cov_matrix = read(fid["covariance_matrix"]) ./ 64000
    close(fid)
    projected_vectors = downsampled_matrix[:, sample_bitmask] * template_projection
    noise_distribution = noise_distribution_given_covaraince_matrix(cov_matrix)
    cur_dist = get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[250, :])
    sum(cur_dist .< cur_dist[intermediate_value_vector[250]])
    plot(get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[250, :]))
    for j in eachindex(intermediate_value_vector)
        most_likely_value = findmax(get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[j, :]))[2] - 1
        if intermediate_value_vector[j] == most_likely_value
            all_correct_counts[template_number] += 1
        end
    end
    all_correct_counts
# end