using HDF5, Base.Threads, StatsBase, Statistics
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")


fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_50_traces_validation.hdf5", "r")

all_intermediate_values = read(fid["intermediate_values"])
downsampled_matrix = read(fid["downsampled_matrix"])
number_of_templates = 2672

for dilation_amount in 0:2:16

    all_correct_counts = zeros(Int64, number_of_templates)
    Threads.@threads for template_number in 1:number_of_templates
        # println(template_number)
        template_path = string("D:\\ChaChaData\\attack_profiling\\initial_templates_infront_50_", dilation_amount, "\\", template_number, "_template.hdf5")
        intermediate_value_vector = all_intermediate_values[:, template_number]
        fid = h5open(template_path, "r")
        sample_bitmask = read(fid["downsampled_sample_bitmask"])
        template_projection = read(fid["projection"])
        mean_vectors = read(fid["class_means"])
        cov_matrix = read(fid["covariance_matrix"])
        close(fid)
        projected_vectors = downsampled_matrix[:, sample_bitmask] * template_projection
        noise_distribution = noise_distribution_given_covaraince_matrix(cov_matrix)
        # val = 12
        # cur_dist = get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[val, :])
        # println(sum(cur_dist .< cur_dist[intermediate_value_vector[val]]))
        # plot(get_prob_dist_of_vector(mean_vectors', noise_distribution, projected_vectors[val, :]))
        for j in eachindex(intermediate_value_vector)
            most_likely_value = findmax(get_prob_dist_of_vector(mean_vectors, noise_distribution, projected_vectors[j, :]))[2] - 1
            if intermediate_value_vector[j] == most_likely_value
                all_correct_counts[template_number] += 1
            end
        end
    end

    # Need to create the output of the success rates as probabilties
    success_rate_percentages = all_correct_counts ./ 10
    # average_percentage_success_rates = success_rate_percentages[begin + 3:4:end]
    average_percentage_success_rates = round.(collect(Iterators.map(mean, Iterators.partition(success_rate_percentages, 4))), digits=2)
    println(dilation_amount, " ", mean(average_percentage_success_rates))
end

# key_success_rates = average_percentage_success_rates[1:8]
# intermediate_value_success_rates = average_percentage_success_rates[29:end]

# function make_quarter_cycle_values(a, b, c, d)
#     return repeat([string(a, "_add"), string(d, "_rot"), string(c, "_add"), string(b, "_rot")], 2)
# end

# function make_break_down_of_values()
#     intermediate_value_locations::Vector{String} = []
#     append!(intermediate_value_locations, make_quarter_cycle_values(1, 5, 9, 13))
#     append!(intermediate_value_locations, make_quarter_cycle_values(2, 6, 10, 14))
#     append!(intermediate_value_locations, make_quarter_cycle_values(3, 7, 11, 15))
#     append!(intermediate_value_locations, make_quarter_cycle_values(4, 8, 12, 16))

#     append!(intermediate_value_locations, make_quarter_cycle_values(1, 6, 11, 16))
#     append!(intermediate_value_locations, make_quarter_cycle_values(2, 7, 12, 13))
#     append!(intermediate_value_locations, make_quarter_cycle_values(3, 8, 9, 14))
#     append!(intermediate_value_locations, make_quarter_cycle_values(4, 5, 10, 15))
#     intermediate_value_locations = repeat(intermediate_value_locations, 10)
#     # append!(intermediate_value_locations, [string(i, "_add") for i in 1:16])
#     return intermediate_value_locations
# end

# intermediate_locations = make_break_down_of_values()

# # Need to break down the cycles into different ways of having counts
# open("average_success_rates_50_4.csv", "w") do file
#     # First do the key
#     write(file, "Key\n")
#     for i in 5:12
#         write(file, string(i, ", ", key_success_rates[i-4], "\n"))
#     end
#     # First do the adds
#     write(file, "Add\n")
#     for i in 1:16
#         if sum(intermediate_locations .== string(i, "_add")) > 0
#             write(file, string(i, ", ", string(intermediate_value_success_rates[intermediate_locations.==string(i, "_add")])[2:end-1], "\n"))
#         end
#     end
#     write(file, "Rotate\n")
#     for i in 1:16
#         if sum(intermediate_locations .== string(i, "_rot")) > 0
#             write(file, string(i, ", ", string(intermediate_value_success_rates[intermediate_locations.==string(i, "_rot")])[2:end-1], "\n"))
#         end
#     end
# end
