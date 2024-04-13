using HDF5, MultivariateStats, Plots, StatsBase, Statistics, LinearAlgebra
include("../encryption/leakage_functions.jl")

# Things which I could try would be to have seperate templats based on the min, mean and max
# to see if they can improve the overall success rate of the template compared to just having a single
# template containing all of them, could also try 'training'

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

# function find_success_rate(distribution, all_positions, mean_vectors, all_values)
#     success_rate = 0
#     @inbounds for i in eachindex(all_values)
#         if all_values[i] == findmax(logpdf(distribution, mean_vectors' .- all_positions[:, i]))[2] - 1
#             success_rate += 1
#         end
#     end
#     return success_rate / length(all_values)
# end

# function calculate_class_means(data_samples, labels)
#     output_matrix = zeros(size(data_samples)[1], length(unique(labels)))
#     for i in unique(labels)
#         output_matrix[:, i + 1] = mean(data_samples[:, labels .== i], dims=2)[:, 1]
#     end
#     return output_matrix
# end

# function calculate_within_class_scatter(data_samples, labels)
#     # Data samples has the dimensions of (d, n) where n is the number of samples and d is the dimensions of the data
#     output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
#     for i in unique(labels)
#         output_matrix += StatsBase.cov(data_samples[labels .== i, :], corrected=true)
#     end
#     return output_matrix ./ length(unique(labels))
# end

# function calculate_between_class_scatter(data_samples, labels, means)
#     # Data samples has the dimensions of (d, n) where n is the number of samples and d is the dimensions of the data
#     output_matrix = zeros(size(data_samples)[1], size(data_samples)[1])
#     overall_mean = mean(matrix_of_current_data, dims=2)[:, 1]
#     for i in unique(labels)
#         vals_to_include = labels .== i
#         cur_mean = means[:, i + 1]
#         output_matrix += ((cur_mean - overall_mean) * (cur_mean - overall_mean)') .* sum(vals_to_include)
#     end
#     return output_matrix ./ length(unique(labels))
# end

# # regularize a symmetric matrix
# function regularize_symmat!(A::AbstractMatrix{T}, lambda::Real) where T<:Real
#     if lambda > 0
#         emax = eigmax(Symmetric(A))
#         add_diag!(A, emax * lambda)
#     end
#     return A
# end

# function add_diag!(A::AbstractMatrix, v::Real)
#     # add v to diagonal of A
#     m = size(A, 1)
#     n = size(A, 2)
#     @assert m == n
#     if v != zero(v)
#         for i = 1:n
#             @inbounds A[i,i] += v
#         end
#     end
#     return A
# end

gr()
number_of_intermediate_values = 2672#32 #2672
number_of_samples_per_cycle = 10
final_number_of_samples = 14989

number_of_output_dimensions = 8
number_of_bits_per_template = 8
number_of_templates_per_intermediate_byte = 8 ÷ number_of_bits_per_template

bitmask_path = "D:\\ChaChaData\\attack_profiling\\clock_cycles_bitmasks_no_dilation.hdf5"
data_path = "D:\\ChaChaData\\attack_profiling\\downsampled_50_traces_profiling.hdf5"
data_path_min = "D:\\ChaChaData\\attack_profiling\\downsampled_50_traces_minimum_profiling.hdf5"
data_path_max = "D:\\ChaChaData\\attack_profiling\\downsampled_50_traces_maximum_profiling.hdf5"
path_to_templates = "D:\\ChaChaData\\attack_profiling\\initial_templates_infront_4_behind_2_min_max_50_"
# intermediate_value_index = 1

# function calculate_within_class_scatter(data_samples, labels, means)
#     # Data samples has the dimensions of (d, n) where n is the number of samples and d is the dimensions of the data
#     output_matrix = zeros(size(data_samples)[1], size(data_samples)[1])
#     for i in eachindex(labels)
#         cur_vec = data_samples[:, i] .- means[labels[i] + 1, :]
#         output_matrix += cur_vec * cur_vec'
#     end
#     return output_matrix ./ length(labels)
# end
# intermediate_value_index = 2
# intermediate_value_index = 1001
# dilation_amount = 0
for dilation_amount in 4:4#0:2:16
    Threads.@threads for intermediate_value_index in 1:number_of_intermediate_values 
        println(intermediate_value_index)
        bitmask_fid = h5open(bitmask_path, "r")
        cycle_bitmask = dilate_after(read(bitmask_fid[string("bitmask_", intermediate_value_index)]), 2)
        cycle_bitmask = dilate_infront(cycle_bitmask, dilation_amount)
        # cycle_bitmask = read(bitmask_fid[string("bitmask_", intermediate_value_index)])
        close(bitmask_fid)

        sample_bitmask = repeat(cycle_bitmask, inner=number_of_samples_per_cycle)[1:final_number_of_samples]
        data_fid = h5open(data_path, "r")
        intermediate_value_vector = read(data_fid["intermediate_values"])[:, intermediate_value_index]
        dset = data_fid["downsampled_matrix"]
        if HDF5.ismmappable(dset)
            dset = HDF5.readmmap(dset)
        end

        # fid = h5open("D:\\ChaChaData\\attack_profiling\\original_1002_values.hdf5", "r")
        # intermediate_value_vector = read(fid["intermediate_values"])
        # original_matrix_of_current_data = fid["downsampled_matrix"][:, 1:3100]'
        # close(fid)


        original_matrix_of_current_data_mean = dset[:, sample_bitmask]'
        close(data_fid)
        
        data_fid = h5open(data_path_min, "r")
        dset = data_fid["downsampled_matrix"]
        if HDF5.ismmappable(dset)
            dset = HDF5.readmmap(dset)
        end
        original_matrix_of_current_data_min = dset[:, sample_bitmask]'
        close(data_fid)
        
        data_fid = h5open(data_path_max, "r")
        dset = data_fid["downsampled_matrix"]
        if HDF5.ismmappable(dset)
            dset = HDF5.readmmap(dset)
        end
        original_matrix_of_current_data_max = dset[:, sample_bitmask]'
        close(data_fid)
        original_matrix_of_current_data = vcat(original_matrix_of_current_data_mean, original_matrix_of_current_data_min, original_matrix_of_current_data_max)
        # println(size(original_matrix_of_current_data))
        template_number = 1
        # for template_number in 1:number_of_templates_per_intermediate_byte
        # println(intermediate_value_index, " ", template_number)
        # Need to get the intermediate value vector for that part of the value
        template_intermediate_value_vector = (intermediate_value_vector .>> (number_of_bits_per_template * (template_number - 1))) .& ((1 << number_of_bits_per_template) - 1)
        permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
        template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
        matrix_of_current_data = original_matrix_of_current_data[:, permutation_of_intermediate_values]

        # original_mean_vectors = calculate_class_means(matrix_of_current_data, template_intermediate_value_vector)
        # test = Matrix(matrix_of_current_data')
        # within_class_scatter = calculate_within_class_scatter(test, template_intermediate_value_vector)
        # between_class_scatter = calculate_between_class_scatter(matrix_of_current_data, template_intermediate_value_vector, original_mean_vectors)

        # within_class_scatter = regularize_symmat!(within_class_scatter, 1e-6)
        # E = eigen!(Symmetric(between_class_scatter), Symmetric(within_class_scatter))
        # ord = sortperm(E.values; rev=true)
        # P = E.vectors[:, ord[1:size(between_class_scatter, 1)]]

        if size(matrix_of_current_data)[1] >= 3000
            lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector; method=:whiten)
        else
            lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector)
        end

        # pca = fit(PCA, matrix_of_current_data)
        # println("test")

        projection_to_subspace = projection(lda)[:, 1:number_of_output_dimensions]
        original_class_means = classmeans(lda)
        within_class_scatter = MultivariateStats.withclass_scatter(lda)
        between_class_scatter = MultivariateStats.betweenclass_scatter(lda)

        projected_class_means = (original_class_means' * projection_to_subspace)
        projected_within_class_scatter = (within_class_scatter * projection_to_subspace)' * projection_to_subspace

        original_data_projected = projection_to_subspace' * matrix_of_current_data
        # projected_scatter = calculate_within_class_scatter(original_data_projected, intermediate_value_vector, projected_class_means)

        # Think my covaraince matrix may be quite wrong really it seems to not have the correct, it is just giving out far to large distriubtion between the 
        # the values which is then resulting in a very weird distribution where pretty much everything is uniform when it really should not
        # be, not really sure where the problem is because I do actually need to get out the projected covariance matrix to be able to model the output well
        # because without it then we are just guessing pretty much at random, like randomly generating values gives much much wider distribution than 
        # what is observed in anyway

        # I believe the covariance matrix needs to be divided by the number of samples which made it up which in this case is 64,000, this seems to give very reasonable 
        # results when actually looking at the scatter in the first couple of dimensions and means we should not see entirely uniform results out when actually doing the templates
        # val = 255
        # current_value_matrix = (matrix_of_current_data[:, template_intermediate_value_vector.==val]' * projection_to_subspace)
        # current_mean = projected_class_means[val+1, :]
        # dist = noise_distribution_given_covaraince_matrix(projected_within_class_scatter)
        # p = scatter(size=(1000, 1000))
        # # randomly_generated_values_around_mean = rand(dist, 1_000) .+ (current_mean)
        # dim_1 = 1
        # dim_2 = 2
        # # scatter!(p, [original_data_projected[1, :]], [original_data_projected[2, :]])
        # scatter!(p, [current_value_matrix[:, dim_1] .* sqrt(64_000)], [current_value_matrix[:, dim_2] .* sqrt(64_000)], label="Actual samples of particular value")
        # # scatter!(p, [randomly_generated_values_around_mean[1, :]], [randomly_generated_values_around_mean[2, :]], label="Generated samples of particular value")
        # scatter!(p, [projected_class_means[:, dim_1] .* sqrt(64_000)], [projected_class_means[:, dim_2] .* sqrt(64_000)], label="Class Means")
        # scatter!(p, [current_mean[dim_1] * sqrt(64_000)], [current_mean[dim_2] * sqrt(64_000)], markersize=10, label="Current class mean")

        # noise_distribution_given_covaraince_matrix()
        if !isdir(string(path_to_templates, dilation_amount))
            mkpath(string(path_to_templates, dilation_amount))
        end
        fid = h5open(string(path_to_templates, dilation_amount, "\\", intermediate_value_index, "_template.hdf5"), "w")
        fid["projection"] = projection_to_subspace .* sqrt(64_000)
        fid["class_means"] = projected_class_means .* sqrt(64_000)
        fid["covariance_matrix"] = projected_within_class_scatter
        fid["downsampled_sample_bitmask"] = sample_bitmask
        close(fid)
    end
end
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