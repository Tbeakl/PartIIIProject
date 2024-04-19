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

function calculate_within_class_scatter(data_samples, labels)
    # Data samples has the dimensions of (n, d) where n is the number of samples and d is the dimensions of the data
    output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
    for i in unique(labels)
        output_matrix += StatsBase.cov(data_samples[labels.==i, :], corrected=true)
    end
    return output_matrix ./ length(unique(labels))
end

function calculate_between_class_scatter(data_samples, labels, means)
    # Data samples has the dimensions of (d, n) where n is the number of samples and d is the dimensions of the data
    output_matrix = zeros(size(data_samples)[1], size(data_samples)[1])
    overall_mean = mean(matrix_of_current_data, dims=2)[:, 1]
    for i in unique(labels)
        vals_to_include = labels .== i
        cur_mean = means[:, i+1]
        output_matrix += ((cur_mean - overall_mean) * (cur_mean - overall_mean)') .* sum(vals_to_include)
    end
    return output_matrix ./ length(unique(labels))
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

gr()
number_of_intermediate_values = 2672 #32 #2672
number_of_samples_per_cycle = 10

number_of_output_dimensions = 8
number_of_bits_per_template = 8
number_of_templates_per_intermediate_byte = 8 ÷ number_of_bits_per_template

bitmask_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation.hdf5"
data_path_20 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_20_traces_profiling.hdf5"
final_number_of_samples_20 = 37471
data_path_50 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_50_traces_profiling.hdf5"
final_number_of_samples_50 = 14989
data_path_100 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_100_traces_profiling.hdf5"
final_number_of_samples_100 = 7495
data_path_250 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_250_traces_profiling.hdf5"
final_number_of_samples_250 = 2998
data_path_500 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_500_traces_profiling.hdf5"
final_number_of_samples_500 = 1499
# data_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_50_traces_maximum_profiling.hdf5"
path_to_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/all_second_templates_LR/"
# intermediate_value_index = 1

function calculate_within_class_scatter(data_samples, labels, means)
    # Data samples has the dimensions of (n, d) where n is the number of samples and d is the dimensions of the data
    output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
    for i in unique(labels)
        part_of_interest = data_samples[labels.==i, :] .- means[i+1, :]'
        output_matrix += Statistics.cov(StatsBase.SimpleCovariance(; corrected=true), part_of_interest; mean=zeros(size(part_of_interest)))
    end
    return output_matrix ./ length(unique(labels))
end

# function calculate_within_class_scatter(data_samples, labels, means)
#     # Data samples has the dimensions of (n, d) where n is the number of samples and d is the dimensions of the data
#     output_matrix = zeros(size(data_samples)[2], size(data_samples)[2])
#     for i in eachindex(labels)
#         println(i)
#         cur_vec = data_samples[i, :] .- means[labels[i]+1, :]
#         output_matrix += cur_vec * cur_vec'
#     end
#     return output_matrix ./ length(labels)
# end
# intermediate_value_index = 2
# intermediate_value_index = 1001
# dilation_amount = 0

data_20_fid = h5open(data_path_20, "r")

dset_intermediate_values = data_20_fid["intermediate_values"]
if HDF5.ismmappable(dset_intermediate_values)
    dset_intermediate_values = HDF5.readmmap(dset_intermediate_values)
end

dset_20 = data_20_fid["downsampled_matrix"]
if HDF5.ismmappable(dset_20)
    dset_20 = HDF5.readmmap(dset_20)
end

data_50_fid = h5open(data_path_50, "r")

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

different_detail_levels = [50, 100, 250, 500]#[20, 50, 100, 250, 500] # [20,50,100,250,500]
different_detail_number_of_samples = [final_number_of_samples_50, final_number_of_samples_100, final_number_of_samples_250, final_number_of_samples_500] #[final_number_of_samples_20, final_number_of_samples_50, final_number_of_samples_100, final_number_of_samples_250, final_number_of_samples_500]
datasets = [dset_50, dset_100, dset_250, dset_500] #[dset_20, dset_50, dset_100, dset_250, dset_500]

simple_intermediate_values = permutedims(hcat(digits.(collect(0:255), base=2, pad=9)...))
simple_intermediate_values[:, end] .= 1
for detailed_level_index in eachindex(different_detail_levels)
    for sparse_level_index in detailed_level_index:length(different_detail_levels)
        println("Detailed: ", detailed_level_index, " Sparse: ", sparse_level_index)
        for intermediate_value_index in 1:number_of_intermediate_values  # Threads.@threads
            if !ispath(string(path_to_templates, "sparse_", different_detail_levels[sparse_level_index], "_detailed_", different_detail_levels[detailed_level_index], "/", intermediate_value_index, "_template.hdf5"))
                println(intermediate_value_index)
                bitmask_fid = h5open(bitmask_path, "r")
                undilated_cycle_bitmask = read(bitmask_fid[string("bitmask_", intermediate_value_index)])
                cycle_bitmask = dilate_after(undilated_cycle_bitmask, 3)
                cycle_bitmask = dilate_infront(cycle_bitmask, 4)
                # cycle_bitmask = read(bitmask_fid[string("bitmask_", intermediate_value_index)])
                close(bitmask_fid)
                intermediate_value_vector = dset_intermediate_values[:, intermediate_value_index]

                detailed_bitmask = repeat(undilated_cycle_bitmask, inner=500 ÷ different_detail_levels[detailed_level_index])[1:different_detail_number_of_samples[detailed_level_index]]
                sparse_bitmask = repeat(cycle_bitmask .⊻ undilated_cycle_bitmask, inner=500 ÷ different_detail_levels[sparse_level_index])[1:different_detail_number_of_samples[sparse_level_index]]
                original_matrix_of_current_data = hcat(datasets[detailed_level_index][:, detailed_bitmask], datasets[sparse_level_index][:, sparse_bitmask])'

                template_number = 1
                template_intermediate_value_vector = (intermediate_value_vector .>> (number_of_bits_per_template * (template_number - 1))) .& ((1 << number_of_bits_per_template) - 1)
                permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
                template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
                matrix_of_current_data = original_matrix_of_current_data[:, permutation_of_intermediate_values]

                # Try and make it using a linear regression model on bits for the mean vectors etc.
                original_mean_vectors = zeros(256, size(matrix_of_current_data)[1])
                base_prediction_matrix = permutedims(hcat(digits.(template_intermediate_value_vector, base=2, pad=9)...))
                base_prediction_matrix[:, end] .= 1
                for i in 1:size(matrix_of_current_data)[1]
                    β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * matrix_of_current_data[i, :]
                    original_mean_vectors[:, i] = simple_intermediate_values * β
                end

                test = Matrix(matrix_of_current_data')

                # Just try a slightly different method for estimating the covariance 

                # Just trying a slightly different method for finding the mean vectors
                # for i in 0:255
                #     original_mean_vectors[i + 1, :] = mean(matrix_of_current_data[:, template_intermediate_value_vector .== i], dims=2)[:, 1]
                # end

                within_class_scatter = calculate_within_class_scatter(test, template_intermediate_value_vector, original_mean_vectors)
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

                projection_to_subspace = P[:, 1:8]
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
                # val = 100
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
                if !isdir(string(path_to_templates, "sparse_", different_detail_levels[sparse_level_index], "_detailed_", different_detail_levels[detailed_level_index], "/"))
                    mkpath(string(path_to_templates, "sparse_", different_detail_levels[sparse_level_index], "_detailed_", different_detail_levels[detailed_level_index], "/"))
                end
                fid = h5open(string(path_to_templates, "sparse_", different_detail_levels[sparse_level_index], "_detailed_", different_detail_levels[detailed_level_index], "/", intermediate_value_index, "_template.hdf5"), "w")
                fid["projection"] = projection_to_subspace
                fid["class_means"] = projected_class_means
                fid["covariance_matrix"] = projected_within_class_scatter
                fid["detailed_sample_bitmask"] = collect(detailed_bitmask)
                fid["sparse_sample_bitmask"] = collect(sparse_bitmask)
                close(fid)
            end
        end
    end
end
close(data_20_fid)
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