using Plots, Base.Threads, Random, HDF5

include("../belief_propagation/node.jl")
include("../belief_propagation/messages.jl")
include("../chacha_factor_graph/chacha_factor_graph.jl")
include("../chacha_factor_graph/add_leakage_to_graph.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")
include("../encryption/leakage_functions.jl")
include("../encryption/chacha.jl")
include("../encryption/rank_estimation.jl")
include("../attacks/fragment_template_attacks/template_attack_traces.jl")

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"
base_key_templates = path_to_data * "attack_profiling/32/initial_templates_sixteen_bit_templates/sparse_50_detailed_50/"
base_path_to_attacks = path_to_data * "captures/ChaChaRecordings_2/recording_attack_counter_from_random_"
base_path_to_fragments = path_to_data * "evaluation/16_bit_fragments_1/"


number_of_bits::Int64 = 16
number_of_templates_per_word::Int64 = 32 ÷ number_of_bits
key_number::Int64 = 1
for key_number in 1:1000
    println(key_number)
    key, nonce, counter, downsampled_trace = load_attack_trace_argmin_alignment(base_path_to_attacks, key_number, 1, path_to_data * "attack_profiling/32/mean_trace.hdf5")

    # Need to make the distributions for each of the parts of the key
    probability_tables::Vector{Vector{Float64}} = []

    for key_word in 1:8
        for template_number in 1:number_of_templates_per_word
            template_path = string(base_key_templates, key_word, "_", template_number, "_template.hdf5")
            fid = h5open(template_path, "r")
            projection_matrix = read(fid["projection"])
            mean_vectors = read(fid["class_means"])
            covaraince_matrix = read(fid["covariance_matrix"])
            sample_bitmask = read(fid["sample_bitmask"])
            close(fid)
            noise = noise_distribution_given_covaraince_matrix(covaraince_matrix)
            position = (downsampled_trace[sample_bitmask]'*projection_matrix)[1, :]
            prob_dist_for_template = get_log_likelihoods_dist_of_vector(mean_vectors, noise, position)
            push!(probability_tables, prob_dist_for_template)
        end
    end

    initial_estimated_rank = rank_estimate_key(key, probability_tables, 500, number_of_bits)
    fid = h5open(string(base_path_to_fragments, key_number, ".hdf5"), "w")
    fid["initial_likelihood_tables"] = reduce(vcat, transpose.(probability_tables))
    initial_estimated_rank = max(1, initial_estimated_rank)
    fid["initial_estimated_rank_log2"] = Float64(log2(initial_estimated_rank))
    close(fid)
end
