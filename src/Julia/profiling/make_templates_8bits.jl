using HDF5, MultivariateStats, Plots, StatsBase, Statistics, LinearAlgebra, FFTW, CUDA
include("../encryption/leakage_functions.jl")
include("common_functions.jl")
# Things which I could try would be to have seperate templats based on the min, mean and max
# to see if they can improve the overall success rate of the template compared to just having a single
# template containing all of them, could also try 'training'
gr()
number_of_intermediate_values = 700

number_of_bits_per_template = 8
number_of_templates_per_intermediate_value = 32 ÷ number_of_bits_per_template
path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"
bitmask_path = path_to_data * "attack_profiling/8_on_32/clock_cycles_bitmasks.hdf5"

data_path = path_to_data * "attack_profiling/8_on_32/"
path_to_templates = path_to_data * "attack_profiling/8_on_32/initial_templates/"

bitmask_fid = h5open(bitmask_path, "r")

function calculate_linear_prediction(base_prediction_matrix, values_per_cycle, simple_intermediate_values, i)
    β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * values_per_cycle[i, :]
    predicted_values = simple_intermediate_values * β
    return predicted_values
end

simple_intermediate_values::Matrix{Int32} = permutedims(hcat(digits.(collect(0:((1<<number_of_bits_per_template)-1)), base=2, pad=number_of_bits_per_template + 1)...))
simple_intermediate_values[:, end] .= 1
simple_intermediate_values = (simple_intermediate_values)

number_of_trace_files = 4

all_data_fids = []
all_data_datasets = []
all_intermediate_datasets = []

for i in 1:number_of_trace_files
    data_fid = h5open(string(data_path, "profiling_2_", i, ".hdf5"), "r")
    dset_intermediate_values = data_fid["intermediate_values"]
    if HDF5.ismmappable(dset_intermediate_values)
        dset_intermediate_values = HDF5.readmmap(dset_intermediate_values)
    end

    dset_data_values = data_fid["downsampled_matrix"]
    if HDF5.ismmappable(dset_data_values)
        dset_data_values = HDF5.readmmap(dset_data_values)
    end

    push!(all_data_fids, data_fid)
    push!(all_data_datasets, dset_data_values)
    push!(all_intermediate_datasets, dset_intermediate_values)
end

Threads.@threads for intermediate_value_index in 1:number_of_intermediate_values
    intermediate_value_vector = make_intermediate_value_matrix(all_intermediate_datasets, intermediate_value_index)
    for template_num in 1:number_of_templates_per_intermediate_value
        current_full_template_path = string(path_to_templates,
            intermediate_value_index, "_", template_num, "_template.hdf5")
        bitmask_dataset_name = string("bitmask_", intermediate_value_index, "_", template_num)
        template_intermediate_value_vector = (intermediate_value_vector .>> (number_of_bits_per_template * (template_num - 1))) .& ((1 << number_of_bits_per_template) - 1)
        permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
        template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
        println(intermediate_value_index, " ", template_num)
        cycle_bitmask = read(bitmask_fid[bitmask_dataset_name])
        if !ispath(current_full_template_path) && sum(cycle_bitmask) > 0
            cycle_bitmask = dilate_infront(cycle_bitmask, 4)
            cycle_bitmask = dilate_after(cycle_bitmask, 2)
            sample_bitmask = Bool.(repeat(cycle_bitmask, inner=50 ÷ 2)[1:399925])
            original_matrix_of_current_data = make_dataset_of_values_matrix(all_data_datasets, sample_bitmask)'

            matrix_of_current_data = original_matrix_of_current_data[:, permutation_of_intermediate_values]

            # Try and make it using a linear regression model on bits for the mean vectors etc.
            original_mean_vectors = zeros(1 << number_of_bits_per_template, size(matrix_of_current_data)[1])
            base_prediction_matrix = ones(Int32, length(template_intermediate_value_vector), number_of_bits_per_template + 1)
            for i in 0:(number_of_bits_per_template-1)
                base_prediction_matrix[:, i+1] = 1 .& (template_intermediate_value_vector .>> i)
            end
            matrix_of_current_data_gpu = (matrix_of_current_data) #CuArray(
            original_mean_vectors_vector = calculate_linear_prediction.(Ref((base_prediction_matrix)), #CuArray(
                Ref(matrix_of_current_data),
                Ref(simple_intermediate_values),
                1:(size(matrix_of_current_data)[1]))
            original_mean_vectors = hcat(original_mean_vectors_vector...)
            template_intermediate_value_vector = (template_intermediate_value_vector)
            within_class_scatter = alternative_within_class_scatter(Matrix(matrix_of_current_data), template_intermediate_value_vector, original_mean_vectors)
            between_class_scatter = alternative_calculate_between_class_scatter(template_intermediate_value_vector, original_mean_vectors)

            within_class_scatter = regularize_symmat!(within_class_scatter, 1e-6)
            E = eigen!(Symmetric(Float64.(between_class_scatter)), Symmetric(Float64.(within_class_scatter)))
            # E = eigen(inv(Symmetric(within_class_scatter))*Symmetric(between_class_scatter)) # Appears to do the same except that it has not been scaled 
            # so the noise has a identity covariance matrix in the projected space
            ord = sortperm(E.values; rev=true)
            P = E.vectors[:, ord[1:size(between_class_scatter, 1)]]

            projection_to_subspace = P[:, 1:number_of_bits_per_template]
            original_class_means = original_mean_vectors

            projected_class_means = (original_class_means * projection_to_subspace)

            original_data_projected = projection_to_subspace' * matrix_of_current_data
            projected_within_class_scatter = alternative_within_class_scatter(Matrix(original_data_projected), template_intermediate_value_vector, projected_class_means)
            # projected_within_class_scatter = calculate_within_class_scatter(original_data_projected', template_intermediate_value_vector, projected_class_means) #(within_class_scatter * projection_to_subspace)' * projection_to_subspace
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
            # dist = noise_distribution_fixed_standard_dev(1., 8) #noise_distribution_given_covaraince_matrix(projected_within_class_scatter)
            # p = scatter(size=(1000, 1000))
            # randomly_generated_values_around_mean = rand(dist, 1_000) .+ (current_mean)
            # dim_1 = 1
            # dim_2 = 2
            # # scatter!(p, [original_data_projected[1, :]], [original_data_projected[2, :]])
            # scatter!(p, [current_value_matrix[:, dim_1]], [current_value_matrix[:, dim_2]], label="Actual samples of particular value")
            # scatter!(p, [randomly_generated_values_around_mean[dim_1, :]], [randomly_generated_values_around_mean[dim_2, :]], label="Generated samples of particular value")
            # scatter!(p, [projected_class_means[:, dim_1]], [projected_class_means[:, dim_2]], label="Class Means")
            # scatter!(p, [current_mean[dim_1]], [current_mean[dim_2]], markersize=10, label="Current class mean")

            # noise_distribution_given_covaraince_matrix()
            if !isdir(string(path_to_templates))
                mkpath(string(path_to_templates))
            end
            fid = h5open(current_full_template_path, "w")
            fid["projection"] = projection_to_subspace
            fid["class_means"] = projected_class_means
            fid["covariance_matrix"] = projected_within_class_scatter
            fid["sample_bitmask"] = collect(sample_bitmask)
            close(fid)
        end
    end
end

for data_fid in all_data_fids
    close(data_fid)
end

close(bitmask_fid)