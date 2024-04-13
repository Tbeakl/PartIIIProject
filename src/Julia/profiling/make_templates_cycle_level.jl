using HDF5, MultivariateStats, Plots, Images
include("../encryption/leakage_functions.jl")

function find_success_rate(distribution, all_positions, mean_vectors, all_values)
    success_rate = 0
    @inbounds for i in eachindex(all_values)
        if all_values[i] == findmax(logpdf(distribution, mean_vectors' .- all_positions[:, i]))[2] - 1
            success_rate += 1
        end
    end
    return success_rate / length(all_values)
end

gr()
number_of_intermediate_values = 1008#32 #2672
number_of_samples_per_cycle = 50
final_number_of_samples = 74941

number_of_output_dimensions = 8
number_of_bits_per_template = 8
number_of_templates_per_intermediate_byte = 8 ÷ number_of_bits_per_template

bitmask_path = "D:\\ChaChaData\\attack_profiling\\clock_cycles_bitmasks_no_dilation.hdf5"
data_path = "D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_profiling.hdf5"
path_to_templates = "D:\\ChaChaData\\attack_profiling\\initial_templates_cycle_level\\"

clock_cycle_sample_number = 405
number_of_intermediate_values = 2800
number_of_clock_cycles = (749500 ÷ 500)

intermediate_values_base_path = "D:\\ChaChaData\\intermediate_value_traces\\recording_profiling_"
traces_base_path = "D:\\ChaChaData\\captures\\ChaChaRecordings\\recording_profiling_"

data_fid = h5open("D:\\ChaChaData\\attack_profiling\\cycle_values_for_NICV.hdf5", "r")
all_intermediate_values = read(data_fid["intermediate_values"])
mean_value_per_cycle = read(data_fid["mean_value_per_cycle"])
min_value_per_cycle = read(data_fid["min_value_per_cycle"])
max_value_per_cycle = read(data_fid["max_value_per_cycle"])
close(data_fid)

data_fid = h5open("D:\\ChaChaData\\attack_profiling\\cycle_values_validation.hdf5", "r")
validation_all_intermediate_values = read(data_fid["intermediate_values"])
validation_mean_value_per_cycle = read(data_fid["mean_value_per_cycle"])
validation_min_value_per_cycle = read(data_fid["min_value_per_cycle"])
validation_max_value_per_cycle = read(data_fid["max_value_per_cycle"])
close(data_fid)

intermediate_value_index = 1933
# for intermediate_value_index in 1001:number_of_intermediate_values
println(intermediate_value_index)
bitmask_fid = h5open(bitmask_path, "r")
cycle_bitmask = read(bitmask_fid[string("bitmask_", intermediate_value_index)])
close(bitmask_fid)

intermediate_value_vector = all_intermediate_values[:, intermediate_value_index]

for dilation_amount in 0:30
    cycle_bitmask_dilated = dilate(cycle_bitmask; r=dilation_amount)
    original_matrix_of_current_data = hcat(mean_value_per_cycle[:, cycle_bitmask_dilated], min_value_per_cycle[:, cycle_bitmask_dilated], max_value_per_cycle[:, cycle_bitmask_dilated])'
    # Need to get the intermediate value vector for that part of the value
    template_intermediate_value_vector = intermediate_value_vector #(intermediate_value_vector .>> (number_of_bits_per_template * (template_number - 1))) .& ((1 << number_of_bits_per_template) - 1)
    permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
    template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
    matrix_of_current_data = original_matrix_of_current_data[:, permutation_of_intermediate_values]

    lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector)

    projection_to_subspace = projection(lda)[:, 1:number_of_output_dimensions]
    original_class_means = classmeans(lda)
    within_class_scatter = MultivariateStats.withclass_scatter(lda)
    between_class_scatter = MultivariateStats.betweenclass_scatter(lda)

    projected_class_means = (original_class_means' * projection_to_subspace)
    projected_within_class_scatter = (within_class_scatter * projection_to_subspace)' * projection_to_subspace

    original_data_projected = projection_to_subspace' * matrix_of_current_data

    validation_data_projected = projection_to_subspace' * hcat(validation_mean_value_per_cycle[:, cycle_bitmask_dilated],
        validation_min_value_per_cycle[:, cycle_bitmask_dilated],
        validation_max_value_per_cycle[:, cycle_bitmask_dilated])'

    distribution = noise_distribution_given_covaraince_matrix(projected_within_class_scatter)
    all_positions = validation_data_projected#[:, 1:1000]
    mean_vectors = projected_class_means
    all_values = validation_all_intermediate_values[:, intermediate_value_index] #[1:1000]
    println(dilation_amount, " ", find_success_rate(distribution, all_positions, mean_vectors, all_values))

    distribution = noise_distribution_given_covaraince_matrix(within_class_scatter)
    all_positions = hcat(validation_mean_value_per_cycle[:, cycle_bitmask_dilated],
    validation_min_value_per_cycle[:, cycle_bitmask_dilated],
    validation_max_value_per_cycle[:, cycle_bitmask_dilated])'
    mean_vectors = original_class_means'
    println(dilation_amount, " ", find_success_rate(distribution, all_positions, mean_vectors, all_values))
end
# projected_scatter = calculate_within_class_scatter(original_data_projected, intermediate_value_vector, projected_class_means)

# Think my covaraince matrix may be quite wrong really it seems to not have the correct, it is just giving out far to large distriubtion between the 
# the values which is then resulting in a very weird distribution where pretty much everything is uniform when it really should not
# be, not really sure where the problem is because I do actually need to get out the projected covariance matrix to be able to model the output well
# because without it then we are just guessing pretty much at random, like randomly generating values gives much much wider distribution than 
# what is observed in anyway

# I believe the covariance matrix needs to be divided by the number of samples which made it up which in this case is 64,000, this seems to give very reasonable 
# results when actually looking at the scatter in the first couple of dimensions and means we should not see entirely uniform results out when actually doing the templates
# val = 0
# current_value_matrix = (matrix_of_current_data[:, template_intermediate_value_vector.==val]' * projection_to_subspace) .* sqrt(64_000)
# current_mean = projected_class_means[val+1, :] .* sqrt(64_000)
# dist = noise_distribution_fixed_standard_dev(1.0, 8) #noise_distribution_given_covaraince_matrix(projected_within_class_scatter)
# p = scatter(size=(1000, 1000))
# randomly_generated_values_around_mean = rand(dist, 1_000) .+ (current_mean)
# dim_1 = 1
# dim_2 = 2
# # scatter!(p, [original_data_projected[1, :]], [original_data_projected[2, :]])
# scatter!(p, [current_value_matrix[:, dim_1]], [current_value_matrix[:, dim_2]], label="Actual samples of particular value")
# # scatter!(p, [randomly_generated_values_around_mean[1, :]], [randomly_generated_values_around_mean[2, :]], label="Generated samples of particular value")
# scatter!(p, [projected_class_means[:, dim_1] .* sqrt(64_000)], [projected_class_means[:, dim_2] .* sqrt(64_000)], label="Class Means")
# scatter!(p, [current_mean[dim_1]], [current_mean[dim_2]], markersize=10, label="Current class mean")

# noise_distribution_given_covaraince_matrix()

# fid = h5open(string(path_to_templates, intermediate_value_index, "_", template_number, "_template.hdf5"), "w")
# fid["projection"] = projection_to_subspace .* sqrt(64_000)
# fid["class_means"] = projected_class_means .* sqrt(64_000)
# fid["covariance_matrix"] = projected_within_class_scatter
# fid["downsampled_sample_bitmask"] = sample_bitmask
# close(fid)
#     end
# end


# Success rate on the training data for a particular template need to find the mean value which has the highest probability for that particular sample


# Currently am losing what order the classes are actually in, which is important
# for actually being able to use the result, I think the most sensible method is to
# sort the values based on the intermediate values so that they are in ascending order so
# hopefully the classes remain in ascending order
# for i in 0:255
#     @assert maximum(abs.((mean(matrix_of_current_data[:, intermediate_value_vector.==i], dims=2)'*projection_to_subspace)[1, :] .- projected_class_means[i+1, :])) < 1e-6
# end