using HDF5, Base.Threads, StatsBase, Statistics, Distributions, LinearAlgebra, Plots
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../chacha_factor_graph/heatmap_visualisation_intermediate_values.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")

function get_prob_dist_of_vector(mean_vectors, noise, current_vector)
    likelihood_of_values = logpdf(noise, mean_vectors' .- current_vector)
    return likelihood_of_values# ./ sum(likelihood_of_values)
end

function get_number_of_guesses(mean_vectors, noise, value, projected_vector)
    permutation_of_results = sortperm(logpdf(noise, mean_vectors .- projected_vector); rev=true)
    return findfirst(x -> x == value + 1, permutation_of_results)
end

function noise_distribution_given_covaraince_matrix(scov)
    return MvNormal(Hermitian(scov))
end

function noise_distribution_fixed_standard_dev(standard_deviation, dimensions)
    variance = standard_deviation * standard_deviation
    covariance_matrix = ones(dimensions) .* variance
    return MvNormal(covariance_matrix)
end

function get_success_rate(means, noise, values, projected_vectors)
    number_of_successes = 0
    for i in eachindex(values)
        most_likely_value = findmax(get_prob_dist_of_vector(means, noise, projected_vectors[i, :]))[2] - 1
        if values[i] == most_likely_value
            number_of_successes += 1
        end
    end
    return number_of_successes / length(values)
end

function get_guessing_entropy(means, noise, values, projected_vectors)
    guessing_entropy = sum(get_number_of_guesses.(Ref(means'), Ref(noise), values, eachrow(projected_vectors)))
    return guessing_entropy / length(values)
end

number_of_bits_per_template = 8
number_of_templates_per_intermediate_value = 32 ÷ number_of_bits_per_template
number_of_itermediate_values = 700# - 32

success_rates = zeros(number_of_itermediate_values * number_of_templates_per_intermediate_value)
guessing_entropies = zeros(number_of_itermediate_values * number_of_templates_per_intermediate_value)

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/8_on_32_trace_set/validation_2.hdf5", "r")
all_intermediate_values = read(fid["intermediate_values"])
downsampled_matrix = read(fid["downsampled_matrix"])
close(fid)

intermediate_value_index = 1
template_number = 1
noise_distribution = noise_distribution_fixed_standard_dev(1.0, number_of_bits_per_template)
for intermediate_value_index in 1:number_of_itermediate_values
    for template_number in 1:number_of_templates_per_intermediate_value
        # Need to change these to work correctly for averaged versions of the templates
        println(intermediate_value_index, " ", template_number)
        template_path = string("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/8_on_32_trace_set/initial_templates/", intermediate_value_index, "_", template_number, "_template.hdf5")
        if ispath(template_path)
            intermediate_value_vector = (all_intermediate_values[:, intermediate_value_index] .>> (number_of_bits_per_template * (template_number - 1))) .& ((1 << number_of_bits_per_template) - 1)
            fid = h5open(template_path, "r")
            cov_matrix = HDF5.readmmap(fid["covariance_matrix"])
            # number_of_dimensions = number_of_bits_per_template# min(number_of_bits_per_template, size(cov_matrix)[1])
            template_projection = HDF5.readmmap(fid["projection"])
            mean_vectors = HDF5.readmmap(fid["class_means"])
            sample_bitmask = HDF5.readmmap(fid["sample_bitmask"])
            # detailed_sample_bitmask = read(fid["detailed_sample_bitmask"])
            # sparse_sample_bitmask = read(fid["sparse_sample_bitmask"])
            original_traces = downsampled_matrix[:, sample_bitmask]
            # original_traces = hcat(downsampled_matrix[:, detailed_sample_bitmask], downsampled_matrix[:, sparse_sample_bitmask])
            projected_vectors = original_traces * template_projection
            # success_rates[number_of_templates_per_intermediate_value * (intermediate_value_index - 1) + template_number] = get_guessing_entropy(mean_vectors, noise_distribution, intermediate_value_vector, projected_vectors)
            success_rates[number_of_templates_per_intermediate_value*(intermediate_value_index-1)+template_number] = get_success_rate(mean_vectors, noise_distribution, intermediate_value_vector, projected_vectors)
            guessing_entropies[number_of_templates_per_intermediate_value*(intermediate_value_index-1)+template_number] = get_guessing_entropy(mean_vectors, noise_distribution, intermediate_value_vector, projected_vectors)
            close(fid)
        end
    end
end

mapping = turn_intermediate_name_to_intermediate_index(number_of_bits_per_template)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits_per_template, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
logarithic_guessing_entopies = map_to_values.(mapped_matrix, Ref(log2.(guessing_entropies)), Ref(2))
# logarithic_guessing_entopies = map_to_values.(mapped_matrix, Ref((success_rates)), Ref(2))

p = heatmap(logarithic_guessing_entopies, title="First order success rate of 16 bit templates", ylabel="Location in state by 16 bit fragments", dpi=300)
# savefig(p, "./plots/heatmaps/16bit_templates_FOSR.png")
# average_percentage_success_rates = round.(collect(Iterators.map(mean, Iterators.partition(success_rates, number_of_templates_per_intermediate_value))), digits=2)
# key_success_rates = average_percentage_success_rates[1:8]
# intermediate_value_success_rates = average_percentage_success_rates[29:end] #-32

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
# open("average_guessing_entropy_constant_attack_sixteen.csv", "w") do file
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

# ge_fid = h5open("guessing_entropies_sixteen_bits_constant_attack_and_success_rates_validation####.hdf5", "w")
# ge_fid["intermediate_success_rate"] = success_rates
# ge_fid["intermediate_guessing_entropies"] = guessing_entropies
# close(ge_fid)