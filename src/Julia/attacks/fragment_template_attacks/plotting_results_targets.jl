using Plots, Base.Threads, Random, HDF5
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("../../encryption/rank_estimation.jl")
include("template_attack_traces.jl")

# initial_key_number = parse(Int64, ARGS[1])
# final_key_number = parse(Int64, ARGS[2])

bits_per_template::Int64 = 8
dimensions_per_template::Int64 = 8
number_of_encryption_traces::Int64 = 10

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

base_path_templates = path_to_data * "attack_profiling/32_volatile/initial_templates_8bits/"
base_key_templates = path_to_data * "attack_profiling/32_volatile/initial_templates_8bits/"
base_trace_path = path_to_data * "captures/ChaChaRecordings_3/recording_attack_counter_constant_"

key_number = 1

all_traces::Matrix{Float32} = zeros(Float32, number_of_encryption_traces, 129960)
key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, 1, path_to_data)
for i in 1:number_of_encryption_traces
    key, nonce, counter, encryption_trace = load_attack_trace(base_trace_path, key_number, i - 1, path_to_data)
    all_traces[i, :] = encryption_trace
end
encryption_output = encrypt(key, nonce, counter)
# There is clearly something not lining up about the intermediate values corresponding to the key or I have really
# badly classified the noise or something, not sure about the other intermediate values because these keys seem way off
elements_of_trace_to_select = append!(repeat([true, false, false, false, true], 320), ones(Bool, 16))
plaintext = zeros(UInt32, 16)
ciphertext = encrypt(key, nonce, counter) .⊻ plaintext
trace = encrypt_collect_trace(key, nonce, counter, full_value)[elements_of_trace_to_select]

all_trace_values = append!(
    key,
    nonce,
    [counter],
    plaintext,
    trace,
    ciphertext)
all_values = UInt32.(collect(Iterators.flatten(all_trace_values)))

intermediate_index_to_plot = 1
template_number_to_use = 1
# profiling_intermediate_value_vector = reduce(vcat, transpose.(byte_values_for_input.(intermediate_values_of_profiling[:, intermediate_index_to_plot])))[:, template_number_to_use]
# traces_to_use_profiling = matrix_used[profiling_intermediate_value_vector.==byte_values_for_input(all_values[intermediate_index_to_plot])[template_number_to_use], :]

p = plot_distribution_of_values_and_means(base_key_templates,
    bits_per_template,
    intermediate_index_to_plot,
    template_number_to_use, byte_values_for_input(all_values[intermediate_index_to_plot])[template_number_to_use],
    all_traces, 1, 2)

mean_trace = median(all_traces[2:10, :], dims=1)[1, :]
p = plot(abs.(all_traces' .- mean_trace), dims=2)
# plot!(p, abs.(all_traces[1, :] .- mean_trace))
plotly()
# savefig(p, "./plots/template_plotting/eight_bits_6_1.svg")