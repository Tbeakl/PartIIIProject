using HDF5, Base.Threads, StatsBase, Statistics
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")

bitmask_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation.hdf5"
data_path_20 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_20_traces_validation.hdf5"
final_number_of_samples_20 = 37471
data_path_50 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_50_traces_validation.hdf5"
final_number_of_samples_50 = 14989
data_path_100 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_100_traces_validation.hdf5"
final_number_of_samples_100 = 7495
data_path_250 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_250_traces_validation.hdf5"
final_number_of_samples_250 = 2998
data_path_500 = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_500_traces_validation.hdf5"
final_number_of_samples_500 = 1499
# data_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_50_traces_maximum_profiling.hdf5"
path_to_templates = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/all_second_templates/"


different_detail_levels = [20, 50, 100, 250, 500] # [20,50,100,250,500]
different_detail_number_of_samples = [final_number_of_samples_20, final_number_of_samples_50, final_number_of_samples_100, final_number_of_samples_250, final_number_of_samples_500]
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

datasets = [dset_20, dset_50, dset_100, dset_250, dset_500]

all_correct_counts = zeros(Int64, length(datasets), length(datasets), number_of_templates)

number_of_templates = 2672

for detailed_level_index in eachindex(different_detail_levels)
    for sparse_level_index in detailed_level_index:length(different_detail_levels)
        println("Detailed: ", detailed_level_index, " Sparse: ", sparse_level_index)
        Threads.@threads for template_number in 1:number_of_templates
            # println(template_number)
            intermediate_value_vector = dset_intermediate_values[:, template_number]
            template_path = string(path_to_templates, "sparse_", different_detail_levels[sparse_level_index], "_detailed_", different_detail_levels[detailed_level_index], "/", template_number, "_template.hdf5")
            fid = h5open(template_path, "r")
            # sample_bitmask = read(fid["downsampled_sample_bitmask"])
            template_projection = read(fid["projection"])
            mean_vectors = read(fid["class_means"])
            cov_matrix = read(fid["covariance_matrix"])
            detailed_sample_bitmask = read(fid["detailed_sample_bitmask"])
            sparse_sample_bitmask = read(fid["sparse_sample_bitmask"])
            close(fid)
            downsampled_matrix = hcat(datasets[detailed_level_index][:, detailed_sample_bitmask], datasets[sparse_level_index][:, sparse_sample_bitmask])
            projected_vectors = downsampled_matrix * template_projection
            noise_distribution = noise_distribution_given_covaraince_matrix(cov_matrix)
            # val = 12
            # cur_dist = get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[val, :])
            # println(sum(cur_dist .< cur_dist[intermediate_value_vector[val]]))
            # plot(get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[val, :]))
            for j in eachindex(intermediate_value_vector)
                most_likely_value = findmax(get_prob_dist_of_vector(mean_vectors, noise_distribution, projected_vectors[j, :]))[2] - 1
                if intermediate_value_vector[j] == most_likely_value
                    all_correct_counts[sparse_level_index, detailed_level_index, template_number] += 1
                end
            end
        end
    end
end

max_performance = maximum(all_correct_counts, dims=[1,2])[1,1,:] ./ 10

# Need to create the output of the success rates as probabilties
success_rate_percentages = all_correct_counts ./ 10
# average_percentage_success_rates = success_rate_percentages[begin+3:4:end]
# average_percentage_success_rates = round.(collect(Iterators.map(mean, Iterators.partition(success_rate_percentages, 4))), digits=2)
average_percentage_success_rates = round.(collect(Iterators.map(mean, Iterators.partition(max_performance, 4))), digits=2)
# println(dilation_amount, ": ", mean(success_rate_percentages))
# println(mean(success_rate_percentages))
# end
key_success_rates = average_percentage_success_rates[1:8]
intermediate_value_success_rates = average_percentage_success_rates[29:end]

function make_quarter_cycle_values(a, b, c, d)
    return repeat([string(a, "_add"), string(d, "_rot"), string(c, "_add"), string(b, "_rot")], 2)
end

function make_break_down_of_values()
    intermediate_value_locations::Vector{String} = []
    append!(intermediate_value_locations, make_quarter_cycle_values(1, 5, 9, 13))
    append!(intermediate_value_locations, make_quarter_cycle_values(2, 6, 10, 14))
    append!(intermediate_value_locations, make_quarter_cycle_values(3, 7, 11, 15))
    append!(intermediate_value_locations, make_quarter_cycle_values(4, 8, 12, 16))

    append!(intermediate_value_locations, make_quarter_cycle_values(1, 6, 11, 16))
    append!(intermediate_value_locations, make_quarter_cycle_values(2, 7, 12, 13))
    append!(intermediate_value_locations, make_quarter_cycle_values(3, 8, 9, 14))
    append!(intermediate_value_locations, make_quarter_cycle_values(4, 5, 10, 15))
    intermediate_value_locations = repeat(intermediate_value_locations, 10)
    # append!(intermediate_value_locations, [string(i, "_add") for i in 1:16])
    return intermediate_value_locations
end

intermediate_locations = make_break_down_of_values()

# Need to break down the cycles into different ways of having counts
open("average_success_rates_second_set_best_1.csv", "w") do file
    # First do the key
    write(file, "Key\n")
    for i in 5:12
        write(file, string(i, ", ", key_success_rates[i-4], "\n"))
    end
    # First do the adds
    write(file, "Add\n")
    for i in 1:16
        if sum(intermediate_locations .== string(i, "_add")) > 0
            write(file, string(i, ", ", string(intermediate_value_success_rates[intermediate_locations.==string(i, "_add")])[2:end-1], "\n"))
        end
    end
    write(file, "Rotate\n")
    for i in 1:16
        if sum(intermediate_locations .== string(i, "_rot")) > 0
            write(file, string(i, ", ", string(intermediate_value_success_rates[intermediate_locations.==string(i, "_rot")])[2:end-1], "\n"))
        end
    end
end

locations_of_maxes = argmax(all_correct_counts, dims=[1,2])[1,1,:]

best_templates_locations = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/second_template_set_best_templates.hdf5"
best_loc_fid = h5open(best_templates_locations, "w")
for i in eachindex(locations_of_maxes)
    best_loc_fid[string("sparse_", i)] = Int64(locations_of_maxes[i][1])
    best_loc_fid[string("detailed_", i)] = Int64(locations_of_maxes[i][2])
end
close(best_loc_fid)