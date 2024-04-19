using HDF5, Base.Threads, StatsBase, Statistics
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")

function downsample_to_mean(trace::AbstractVector, number_of_samples::Int64)
    return collect(Iterators.map(mean, Iterators.partition(trace, number_of_samples)))
end

function downsample_to_mean_weighted_by_NICV(trace::AbstractVector, nicv::Vector{Float64}, number_of_samples::Int64)
    return collect(Iterators.map(sum, Iterators.partition(trace .* nicv, number_of_samples))) ./ collect(Iterators.map(sum, Iterators.partition(nicv, number_of_samples)))
end

base_intermediate_value_index = 1001

fid = h5open(string("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/validation_", base_intermediate_value_index, ".hdf5"), "r")
all_intermediate_values = read(fid["intermediate_values"])
original_matrix = read(fid["downsampled_matrix"])
close(fid)

bitmask_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation.hdf5"
bitmask_fid = h5open(bitmask_path, "r")
detailed_bitmask = read(bitmask_fid[string("bitmask_", base_intermediate_value_index)])
full_cycle_bitmask = dilate_after(detailed_bitmask, 3)
full_cycle_bitmask = dilate_infront(full_cycle_bitmask, 4)
detailed_bitmask = repeat(detailed_bitmask[full_cycle_bitmask], inner=500)
close(bitmask_fid)

interesting_samples = original_matrix[:, detailed_bitmask]
dilation_samples = original_matrix[:, .!detailed_bitmask]

success_rates_mean = zeros(Int64, 3, 9)
success_rates_nicv = zeros(Int64, 3, 9)

for intermediate_index in base_intermediate_value_index:(base_intermediate_value_index+3)
    println("####################################")
    println(intermediate_index)
    for (dilation_number, number_of_samples_dilation_averaged_over) in enumerate([100, 250, 500])
        downsampled_dilation = reduce(vcat, transpose.(downsample_to_mean.(eachrow(dilation_samples), number_of_samples_dilation_averaged_over)))
        for (i, number_of_samples_to_mean_over) in enumerate([1, 2, 5, 10, 20, 50, 100, 250, 500]) #1, 2, 5, 
            detailed_downsampled_part = reduce(vcat, transpose.(downsample_to_mean.(eachrow(interesting_samples), number_of_samples_to_mean_over)))
            downsampled_matrix = hcat(downsampled_dilation, detailed_downsampled_part)
            all_correct_counts = 0
            template_path = string("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/initial_templates_all_cycles_", intermediate_index, "/dilated_", number_of_samples_dilation_averaged_over, "_", number_of_samples_to_mean_over, "_points_template.hdf5")
            intermediate_value_vector = all_intermediate_values[:, intermediate_index-base_intermediate_value_index+1]
            fid = h5open(template_path, "r")
            template_projection = read(fid["projection"])
            mean_vectors = read(fid["class_means"])
            cov_matrix = read(fid["covariance_matrix"])
            close(fid)

            projected_vectors = downsampled_matrix * template_projection
            noise_distribution = noise_distribution_given_covaraince_matrix(cov_matrix)
            for j in eachindex(intermediate_value_vector)
                most_likely_value = findmax(get_prob_dist_of_vector(mean_vectors, noise_distribution, projected_vectors[j, :]))[2] - 1
                if intermediate_value_vector[j] == most_likely_value
                    success_rates_mean[dilation_number, i] += 1
                end
            end

            # success_rate_percentages = all_correct_counts / 10
            # println("Mean: ", number_of_samples_to_mean_over, ": ", success_rate_percentages)

            all_correct_counts = 0
            template_path = string("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/initial_templates_all_cycles_", intermediate_index, "/dilated_", number_of_samples_dilation_averaged_over, "_", number_of_samples_to_mean_over, "_weighted_points_template.hdf5")
            intermediate_value_vector = all_intermediate_values[:, intermediate_index-base_intermediate_value_index+1]
            fid = h5open(template_path, "r")
            template_projection = read(fid["projection"])
            mean_vectors = read(fid["class_means"])
            cov_matrix = read(fid["covariance_matrix"])
            nicv = read(fid["NICV"])
            close(fid)
            detailed_downsampled_part = reduce(vcat, transpose.(downsample_to_mean_weighted_by_NICV.(eachrow(interesting_samples), Ref(nicv), number_of_samples_to_mean_over)))
            downsampled_matrix = hcat(downsampled_dilation, detailed_downsampled_part)

            projected_vectors = downsampled_matrix * template_projection
            noise_distribution = noise_distribution_given_covaraince_matrix(cov_matrix)
            for j in eachindex(intermediate_value_vector)
                most_likely_value = findmax(get_prob_dist_of_vector(mean_vectors, noise_distribution, projected_vectors[j, :]))[2] - 1
                if intermediate_value_vector[j] == most_likely_value
                    success_rates_nicv[dilation_number, i] += 1
                end
            end

            # success_rate_percentages = all_correct_counts / 10
            # println("NICV: ", number_of_samples_to_mean_over, ": ", success_rate_percentages)
        end
    end
end

# println("Mean: ", success_rates_mean ./ 40)
# println("NICV: ", success_rates_nicv ./ 40)