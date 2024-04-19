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

function calculate_NICV(traces, variances, intermediate_value_vector)
    all_mean_values = zeros(size(traces)[2], 256)
    for j in 0:255
        all_mean_values[:, j+1] = mean(traces[intermediate_value_vector.==j, :], dims=1)[1, :]
    end
    return var(all_mean_values, dims=2) ./ variances
end

function downsample_to_mean(trace::AbstractVector, number_of_samples::Int64)
    return collect(Iterators.map(mean, Iterators.partition(trace, number_of_samples)))
end

function downsample_to_mean_weighted_by_NICV(trace::AbstractVector, nicv::Vector{Float64}, number_of_samples::Int64)
    return collect(Iterators.map(sum, Iterators.partition(trace .* nicv, number_of_samples))) ./ collect(Iterators.map(sum, Iterators.partition(nicv, number_of_samples)))
end

number_of_intermediate_values = 4
base_intermediate_value_index = 1

number_of_output_dimensions = 8
number_of_bits_per_template = 8

bitmask_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation.hdf5"
data_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_25_traces_profiling.hdf5"
path_to_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/initial_templates_all_cycles_"

data_fid = h5open(data_path, "r")

dset_intermediate_values = data_fid["intermediate_values"]
if HDF5.ismmappable(dset_intermediate_values)
    dset_intermediate_values = HDF5.readmmap(dset_intermediate_values)
end

dset = data_fid["downsampled_matrix"]
if HDF5.ismmappable(dset)
    dset = HDF5.readmmap(dset)
end

# variances = var(dset[:, :], dims=1)[1, :]
# all_NICVS = zeros(length(variances), 4)
# for i in 1:4
#     println(i)
#     all_NICVS[:, i] = calculate_NICV(dset[:, :], variances, dset_intermediate_values[:, i])[:, 1]
# end

bitmask_fid = h5open(bitmask_path, "r")
detailed_bitmask = read(bitmask_fid[string("bitmask_", base_intermediate_value_index)])
full_cycle_bitmask = dilate_after(detailed_bitmask, 3)
full_cycle_bitmask = dilate_infront(full_cycle_bitmask, 4)
detailed_bitmask = repeat(detailed_bitmask[full_cycle_bitmask], inner=500)
close(bitmask_fid)

interesting_samples = dset[:, detailed_bitmask]
dilation_samples = dset[:, .!detailed_bitmask]

# number_of_samples_averaged = 2
intermediate_value_index = 1
for dilation_samples_averaging_samples in collect([100, 250, 500])
    downsampled_dilation = reduce(vcat, transpose.(downsample_to_mean.(eachrow(dilation_samples), dilation_samples_averaging_samples)))
    for number_of_samples_averaged in collect([10, 20, 50, 100, 250, 500])
        println(number_of_samples_averaged)
        for intermediate_value_index in 1:number_of_intermediate_values
            intermediate_value_vector = dset_intermediate_values[:, intermediate_value_index]
            # NICVs = all_NICVS[detailed_bitmask, intermediate_value_index]
            # println("Before reducing in shape")
            # detailed_downsampled_part = reduce(vcat, transpose.(downsample_to_mean_weighted_by_NICV.(eachrow(interesting_samples), Ref(NICVs), number_of_samples_averaged)))
            # original_matrix_of_current_data = hcat(downsampled_dilation, detailed_downsampled_part)
            # println(size(original_matrix_of_current_data))
            # if !ispath(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/dilated_", dilation_samples_averaging_samples, "_", number_of_samples_averaged, "_weighted_points_template.hdf5"))
            #     template_intermediate_value_vector = intermediate_value_vector
            #     permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
            #     template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
            #     matrix_of_current_data = Float32.(original_matrix_of_current_data'[:, permutation_of_intermediate_values])

            #     if size(matrix_of_current_data)[1] >= 3000
            #         lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector; method=:whiten)
            #     else
            #         lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector)
            #     end

            #     projection_to_subspace = projection(lda)[:, 1:number_of_output_dimensions]
            #     original_class_means = classmeans(lda)
            #     within_class_scatter = MultivariateStats.withclass_scatter(lda)
            #     between_class_scatter = MultivariateStats.betweenclass_scatter(lda)

            #     projected_class_means = (original_class_means' * projection_to_subspace)
            #     projected_within_class_scatter = (within_class_scatter * projection_to_subspace)' * projection_to_subspace

            #     original_data_projected = projection_to_subspace' * matrix_of_current_data

            #     if !isdir(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/"))
            #         mkpath(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/"))
            #     end
            #     fid = h5open(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/dilated_",dilation_samples_averaging_samples, "_", number_of_samples_averaged, "_weighted_points_template.hdf5"), "w")
            #     fid["projection"] = projection_to_subspace .* sqrt(64_000)
            #     fid["class_means"] = projected_class_means .* sqrt(64_000)
            #     fid["covariance_matrix"] = projected_within_class_scatter
            #     fid["NICV"] = NICVs
            #     # fid["downsampled_sample_bitmask"] = sample_bitmask
            #     close(fid)
            # end

            detailed_downsampled_part = reduce(vcat, transpose.(downsample_to_mean.(eachrow(interesting_samples), number_of_samples_averaged)))
            original_matrix_of_current_data = hcat(downsampled_dilation, detailed_downsampled_part)
            println(size(original_matrix_of_current_data))
            if !ispath(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/dilated_",dilation_samples_averaging_samples, "_", number_of_samples_averaged, "_points_template.hdf5"))

                template_intermediate_value_vector = intermediate_value_vector
                permutation_of_intermediate_values = sortperm(template_intermediate_value_vector)
                template_intermediate_value_vector = template_intermediate_value_vector[permutation_of_intermediate_values]
                matrix_of_current_data = Float32.(original_matrix_of_current_data'[:, permutation_of_intermediate_values])

                if size(matrix_of_current_data)[1] >= 3000
                    lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector; method=:whiten)
                else
                    lda = fit(MulticlassLDA, matrix_of_current_data, template_intermediate_value_vector)
                end

                projection_to_subspace = projection(lda)[:, 1:number_of_output_dimensions]
                original_class_means = classmeans(lda)
                within_class_scatter = MultivariateStats.withclass_scatter(lda)
                between_class_scatter = MultivariateStats.betweenclass_scatter(lda)

                projected_class_means = (original_class_means' * projection_to_subspace)
                projected_within_class_scatter = (within_class_scatter * projection_to_subspace)' * projection_to_subspace

                original_data_projected = projection_to_subspace' * matrix_of_current_data

                if !isdir(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/"))
                    mkpath(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/"))
                end
                fid = h5open(string(path_to_templates, intermediate_value_index + base_intermediate_value_index - 1, "/dilated_",dilation_samples_averaging_samples, "_", number_of_samples_averaged, "_points_template.hdf5"), "w")
                fid["projection"] = projection_to_subspace .* sqrt(64_000)
                fid["class_means"] = projected_class_means .* sqrt(64_000)
                fid["covariance_matrix"] = projected_within_class_scatter
                # fid["NICV"] = NICVs
                # fid["downsampled_sample_bitmask"] = sample_bitmask
                close(fid)
            end
        end
    end
end
close(data_fid)