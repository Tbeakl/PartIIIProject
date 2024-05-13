using HDF5, Base.Threads, StatsBase, Statistics, Distributions, LinearAlgebra, Plots
gr()
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../chacha_factor_graph/heatmap_visualisation_intermediate_values.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")
include("common_functions.jl")

number_of_bits_per_template = 8
number_of_templates_per_intermediate_value = 32 ÷ number_of_bits_per_template
number_of_itermediate_values = 700# - 32

success_rates = zeros(number_of_itermediate_values * number_of_templates_per_intermediate_value)
guessing_entropies = zeros(number_of_itermediate_values * number_of_templates_per_intermediate_value)

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

fid = h5open(path_to_data * "attack_profiling/8_on_32/validation.hdf5", "r")
all_intermediate_values = read(fid["intermediate_values"])
downsampled_matrix = read(fid["downsampled_matrix"])
close(fid)

for intermediate_value_index in 1:number_of_itermediate_values
    for template_number in 1:number_of_templates_per_intermediate_value
        # Need to change these to work correctly for averaged versions of the templates
        println(intermediate_value_index, " ", template_number)
        template_path = string(path_to_data * "attack_profiling/8_on_32/initial_templates/", intermediate_value_index, "_", template_number, "_template.hdf5")
        if ispath(template_path)
            intermediate_value_vector = (all_intermediate_values[:, intermediate_value_index] .>> (number_of_bits_per_template * (template_number - 1))) .& ((1 << number_of_bits_per_template) - 1)
            fid = h5open(template_path, "r")
            cov_matrix = HDF5.readmmap(fid["covariance_matrix"])
            template_projection = HDF5.readmmap(fid["projection"])
            mean_vectors = HDF5.readmmap(fid["class_means"])
            sample_bitmask = HDF5.readmmap(fid["sample_bitmask"])
            original_traces = downsampled_matrix[:, sample_bitmask]
            projected_vectors = original_traces * template_projection
            noise = noise_distribution_given_covaraince_matrix(cov_matrix)
            success_rates[number_of_templates_per_intermediate_value * (intermediate_value_index - 1) + template_number] = get_success_rate(mean_vectors, noise, intermediate_value_vector, projected_vectors)
            guessing_entropies[number_of_templates_per_intermediate_value * (intermediate_value_index - 1) + template_number] = get_guessing_entropy(mean_vectors, noise, intermediate_value_vector, projected_vectors)
            close(fid)
        end
    end
end

mapping = turn_intermediate_name_to_intermediate_index(number_of_bits_per_template)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits_per_template, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
logarithic_guessing_entopies = map_to_values.(mapped_matrix, Ref(log2.(guessing_entropies)), Ref(32 ÷ number_of_bits_per_template))
success_rates_heatmap = map_to_values.(mapped_matrix, Ref(success_rates), Ref(32 ÷ number_of_bits_per_template))
p = heatmap(logarithic_guessing_entopies, title="Logarithmic guessing entropy of byte templates", ylabel="Location in state by bytes", dpi=300)
savefig(p, "./plots/heatmaps/byte_templates_LGE.png")

p = heatmap(success_rates_heatmap, title="First-order success rates of byte templates", ylabel="Location in state by bytes", dpi=300)
savefig(p, "./plots/heatmaps/byte_templates_FSR.png")
