using HDF5, MultivariateStats, Plots, StatsBase, Statistics, LinearAlgebra
include("../encryption/leakage_functions.jl")
include("common_functions.jl")
# Things which I could try would be to have seperate templats based on the min, mean and max
# to see if they can improve the overall success rate of the template compared to just having a single
# template containing all of them, could also try 'training'
gr()
number_of_intermediate_values = 700 #32 #2672

# number_of_output_dimensions = 8
total_variation_to_include = 0.9
number_of_bits_per_template = 8
number_of_templates_per_intermediate_value = 32 ÷ number_of_bits_per_template

bitmask_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_dilated_by_leakage_event.hdf5"

# data_path_20 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/profiling_20.hdf5"
# final_number_of_samples_20 = 37471
data_path_50 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/profiling_50.hdf5"
final_number_of_samples_50 = 14989
data_path_100 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/profiling_100.hdf5"
final_number_of_samples_100 = 7495
data_path_250 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/profiling_250.hdf5"
final_number_of_samples_250 = 2998
data_path_500 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/profiling_500.hdf5"
final_number_of_samples_500 = 1499
# data_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_50_traces_maximum_profiling.hdf5"
path_to_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_trace_set/initial_templates_LR_leakage_event/"

bitmask_fid = h5open(bitmask_path, "r")

# data_20_fid = h5open(data_path_20, "r")

# dset_20 = data_20_fid["downsampled_matrix"]
# if HDF5.ismmappable(dset_20)
#     dset_20 = HDF5.readmmap(dset_20)
# end

data_50_fid = h5open(data_path_50, "r")

dset_intermediate_values = data_50_fid["intermediate_values"]
if HDF5.ismmappable(dset_intermediate_values)
    dset_intermediate_values = HDF5.readmmap(dset_intermediate_values)
end

dset_50 = data_50_fid["downsampled_matrix"]
if HDF5.ismmappable(dset_50)
    dset_50 = HDF5.readmmap(dset_50)
end

data_100_fid = h5open(data_path_100, "r")

dset_100 = data_100_fid["downsampled_matrix"]
if HDF5.ismmappable(dset_100)
    dset_100 = HDF5.readmmap(dset_100)
end

data_250_fid = h5open(data_path_250, "r")

dset_250 = data_250_fid["downsampled_matrix"]
if HDF5.ismmappable(dset_250)
    dset_250 = HDF5.readmmap(dset_250)
end

data_500_fid = h5open(data_path_500, "r")

dset_500 = data_500_fid["downsampled_matrix"]
if HDF5.ismmappable(dset_500)
    dset_500 = HDF5.readmmap(dset_500)
end

different_detail_levels = [50]
different_detail_number_of_samples = [final_number_of_samples_50]
datasets = [dset_50]

simple_intermediate_values = permutedims(hcat(digits.(collect(0:((1<<number_of_bits_per_template)-1)), base=2, pad=number_of_bits_per_template + 1)...))
simple_intermediate_values[:, end] .= 1

# detailed_level_index = 1
# sparse_level_index = 1
# intermediate_value_index = 1
# template_num = 1
# leakage_event_number = 4
detailed_level_index = 1
sparse_level_index = 1
# for detailed_level_index in eachindex(different_detail_levels)
#     for sparse_level_index in detailed_level_index:length(different_detail_levels)
println("Detailed: ", detailed_level_index, " Sparse: ", sparse_level_index)
for intermediate_value_index in 1:number_of_intermediate_values
    intermediate_value_vector = dset_intermediate_values[:, intermediate_value_index]
    for template_num in 1:number_of_templates_per_intermediate_value
        leakage_event_number::Int64 = 1
        current_full_template_path = string(path_to_templates, "sparse_",
            different_detail_levels[sparse_level_index], "_detailed_",
            different_detail_levels[detailed_level_index], "/",
            intermediate_value_index, "_", template_num, "_", leakage_event_number, "_template.hdf5")
        bitmask_dataset_name = string("bitmask_", intermediate_value_index, "_", leakage_event_number)
        template_intermediate_value_vector = (intermediate_value_vector .>> (number_of_bits_per_template * (template_num - 1))) .& ((1 << number_of_bits_per_template) - 1)
        permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
        template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
        while !ispath(current_full_template_path) && haskey(bitmask_fid, bitmask_dataset_name)
            println(intermediate_value_index, " ", template_num, " ", leakage_event_number)
            cycle_bitmask = read(bitmask_fid[bitmask_dataset_name])
            sample_bitmask = Bool.(repeat(cycle_bitmask, inner=500 ÷ different_detail_levels[detailed_level_index])[1:different_detail_number_of_samples[detailed_level_index]])
            original_matrix_of_current_data = datasets[detailed_level_index][:, sample_bitmask]'
            # detailed_bitmask = repeat(undilated_cycle_bitmask, inner=500 ÷ different_detail_levels[detailed_level_index])[1:different_detail_number_of_samples[detailed_level_index]]
            # sparse_bitmask = repeat(cycle_bitmask .⊻ undilated_cycle_bitmask, inner=500 ÷ different_detail_levels[sparse_level_index])[1:different_detail_number_of_samples[sparse_level_index]]
            # original_matrix_of_current_data = hcat(datasets[detailed_level_index][:, detailed_bitmask], datasets[sparse_level_index][:, sparse_bitmask])'

            matrix_of_current_data = original_matrix_of_current_data[:, permutation_of_intermediate_values]

            # Try and make it using a linear regression model on bits for the mean vectors etc.
            original_mean_vectors = zeros(1 << number_of_bits_per_template, size(matrix_of_current_data)[1])
            base_prediction_matrix = ones(length(template_intermediate_value_vector), number_of_bits_per_template + 1)
            for i in 0:(number_of_bits_per_template-1)
                base_prediction_matrix[:, i+1] = 1 .& (template_intermediate_value_vector .>> i)
            end
            for i in 1:size(matrix_of_current_data)[1]
                β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * matrix_of_current_data[i, :]
                original_mean_vectors[:, i] = simple_intermediate_values * β
            end

            transposed_data = Matrix(matrix_of_current_data')

            # Just try a slightly different method for estimating the covariance

            # Just trying a slightly different method for finding the mean vectors
            # for i in 0:255
            #     original_mean_vectors[i + 1, :] = mean(matrix_of_current_data[:, template_intermediate_value_vector .== i], dims=2)[:, 1]
            # end

            # original_mean_vectors = zeros(number_of_bits_per_template + 1, size(matrix_of_current_data)[1])
            # for i in 0:number_of_bits_per_template
            #     original_mean_vectors[i+1, :] = mean(matrix_of_current_data[:, template_intermediate_value_vector.==i], dims=2)[:, 1]
            # end

            within_class_scatter = calculate_within_class_scatter(transposed_data, template_intermediate_value_vector, original_mean_vectors)
            between_class_scatter = calculate_between_class_scatter(matrix_of_current_data, template_intermediate_value_vector, original_mean_vectors')

            within_class_scatter = regularize_symmat!(within_class_scatter, 1e-6)
            E = eigen!(Symmetric(between_class_scatter), Symmetric(within_class_scatter))
            ord = sortperm(E.values; rev=true)
            P = E.vectors[:, ord[1:size(between_class_scatter, 1)]]

            # if size(matrix_of_current_data)[1] >= 3000
            #     lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector; method=:whiten)
            # else
            #     lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector)
            # end

            # pca = fit(PCA, matrix_of_current_data)
            # println("test")

            # # # projection_to_subspace = projection(lda)[:, 1:number_of_output_dimensions]
            # # # original_class_means = classmeans(lda)
            # # # within_class_scatter = MultivariateStats.withclass_scatter(lda)
            # # # between_class_scatter = MultivariateStats.betweenclass_scatter(lda)
            number_of_output_dimensions = number_of_bits_per_template # findfirst(x->x > total_variation_to_include, cumsum(abs.(E.values[ord]) ./ sum(abs.(E.values))))
            projection_to_subspace = P[:, 1:number_of_output_dimensions]
            original_class_means = original_mean_vectors

            projected_class_means = (original_class_means * projection_to_subspace)

            original_data_projected = projection_to_subspace' * matrix_of_current_data
            projected_within_class_scatter = calculate_within_class_scatter(original_data_projected', template_intermediate_value_vector, projected_class_means) #(within_class_scatter * projection_to_subspace)' * projection_to_subspace
            # projected_scatter = calculate_within_class_scatter(original_data_projected, intermediate_value_vector, projected_class_means)

            # Think my covaraince matrix may be quite wrong really it seems to not have the correct, it is just giving out far to large distriubtion between the
            # the values which is then resulting in a very weird distribution where pretty much everything is uniform when it really should not
            # be, not really sure where the problem is because I do actually need to get out the projected covariance matrix to be able to model the output well
            # because without it then we are just guessing pretty much at random, like randomly generating values gives much much wider distribution than
            # what is observed in anyway

            # I believe the covariance matrix needs to be divided by the number of samples which made it up which in this case is 64,000, this seems to give very reasonable
            # results when actually looking at the scatter in the first couple of dimensions and means we should not see entirely uniform results out when actually doing the templates
            # val = 0
            # current_value_matrix = (matrix_of_current_data[:, template_intermediate_value_vector.==val]' * projection_to_subspace)
            # current_mean = projected_class_means[val+1, :]
            # # dist = noise_distribution_given_covaraince_matrix(projected_within_class_scatter)
            # p = scatter(size=(1000, 1000))
            # # randomly_generated_values_around_mean = rand(dist, 1_000) .+ (current_mean)
            # dim_1 = 1
            # dim_2 = 2
            # # scatter!(p, [original_data_projected[1, :]], [original_data_projected[2, :]])
            # scatter!(p, [current_value_matrix[:, dim_1]], [current_value_matrix[:, dim_2]], label="Actual samples of particular value")
            # # scatter!(p, [randomly_generated_values_around_mean[1, :]], [randomly_generated_values_around_mean[2, :]], label="Generated samples of particular value")
            # scatter!(p, [projected_class_means[:, dim_1]], [projected_class_means[:, dim_2]], label="Class Means")
            # scatter!(p, [current_mean[dim_1]], [current_mean[dim_2]], markersize=10, label="Current class mean")

            # noise_distribution_given_covaraince_matrix()
            if !isdir(string(path_to_templates, "sparse_",
                different_detail_levels[sparse_level_index], "_detailed_",
                different_detail_levels[detailed_level_index], "/"))
                mkpath(string(path_to_templates, "sparse_",
                    different_detail_levels[sparse_level_index], "_detailed_",
                    different_detail_levels[detailed_level_index], "/"))
            end
            fid = h5open(current_full_template_path, "w")
            fid["projection"] = projection_to_subspace
            fid["class_means"] = projected_class_means
            fid["covariance_matrix"] = projected_within_class_scatter
            fid["sample_bitmask"] = collect(sample_bitmask)
            close(fid)
            leakage_event_number += 1
            current_full_template_path = string(path_to_templates, "sparse_",
                different_detail_levels[sparse_level_index], "_detailed_",
                different_detail_levels[detailed_level_index], "/",
                intermediate_value_index, "_", template_num, "_", leakage_event_number, "_template.hdf5")
            bitmask_dataset_name = string("bitmask_", intermediate_value_index, "_", leakage_event_number)
        end
    end
end
#     end
# end
# close(data_20_fid)
close(data_50_fid)
close(data_100_fid)
close(data_250_fid)
close(data_500_fid)
# fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_validation.hdf5", "r")


# validation_intermediate_values = read(fid["intermediate_values"])[:, intermediate_value_index]
# dset = fid["downsampled_matrix"]
# if HDF5.ismmappable(dset)
#     dset = HDF5.readmmap(dset)
# end
# downsampled_matrix = dset[:, sample_bitmask]'
# close(fid)
# distribution = noise_distribution_given_covaraince_matrix(projected_within_class_scatter)
# all_positions = (projection_to_subspace' .* sqrt(64_000)) * downsampled_matrix #[:, 1:1000]
# mean_vectors = projected_class_means .* sqrt(64_000)
# all_values = validation_intermediate_values #[1:1000]
# println(dilation_amount, " ", find_success_rate(distribution, all_positions, mean_vectors, all_values))

# distribution = noise_distribution_given_covaraince_matrix(within_class_scatter)
# all_positions = matrix_of_current_data[:, 1:100]
# mean_vectors = original_class_means'
# all_values = template_intermediate_value_vector[1:100]
# find_success_rate(distribution, all_positions, mean_vectors, all_values)
# Success rate on the training data for a particular template need to find the mean value which has the highest probability for that particular sample


# Currently am losing what order the classes are actually in, which is important
# for actually being able to use the result, I think the most sensible method is to
# sort the values based on the intermediate values so that they are in ascending order so
# hopefully the classes remain in ascending order
# for i in 0:255
#     @assert maximum(abs.((mean(matrix_of_current_data[:, intermediate_value_vector.==i], dims=2)'*projection_to_subspace)[1, :] .- projected_class_means[i+1, :])) < 1e-6
# end