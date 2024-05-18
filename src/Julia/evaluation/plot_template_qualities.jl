using Plots, HDF5, Statistics
include("../chacha_factor_graph/heatmap_visualisation_intermediate_values.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")
# 8 bit implementation
fid = h5open("./data/evaluation/heatmap_data/FSR/8_on_32.hdf5", "r")
success_rates = read(fid["success_rates"])
close(fid)
number_of_bits = 8
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(success_rates), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="First-order success rate", ylabel="ChaCha state", xlabel="Operation number", dpi=300)
savefig(p, "./plots/heatmaps/byte_templates_FSR.png")

fid = h5open("./data/evaluation/heatmap_data/LGE/8_on_32.hdf5", "r")
guessing_entropies = read(fid["guessing_entropies"])
close(fid)
number_of_bits = 8
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(log2.(guessing_entropies)), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="Logarithmic guessing entropy", ylabel="ChaCha state", xlabel="Operation number", clim=(0,7), dpi=300)
savefig(p, "./plots/heatmaps/byte_templates_LGE.png")
println("8: LGE: ", log2(mean(guessing_entropies)), " FSR: ", mean(success_rates))

# 32 volatile
fid = h5open("./data/evaluation/heatmap_data/FSR/32_volatile.hdf5", "r")
success_rates = read(fid["success_rates"])
close(fid)
number_of_bits = 8
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(success_rates), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="First-order success rate", ylabel="ChaCha state", xlabel="Operation number", dpi=300)
savefig(p, "./plots/heatmaps/32_volatile_byte_templates_FSR.png")

fid = h5open("./data/evaluation/heatmap_data/LGE/32_volatile.hdf5", "r")
guessing_entropies = read(fid["guessing_entropies"])
close(fid)
number_of_bits = 8
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(log2.(guessing_entropies)), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="Logarithmic guessing entropy", ylabel="ChaCha state", xlabel="Operation number", clim=(0,7), dpi=300)
savefig(p, "./plots/heatmaps/32_volatile_byte_templates_LGE.png")
println("32 volatile: LGE: ", log2(mean(guessing_entropies)), " FSR: ", mean(success_rates))

# 32
fid = h5open("./data/evaluation/heatmap_data/FSR/32.hdf5", "r")
success_rates = read(fid["success_rates"])
close(fid)
number_of_bits = 16
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(success_rates), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="First-order success rate", ylabel="ChaCha state", xlabel="Operation number", dpi=300)
savefig(p, "./plots/heatmaps/32_16_bit_templates_FSR.png")

fid = h5open("./data/evaluation/heatmap_data/LGE/32.hdf5", "r")
guessing_entropies = read(fid["guessing_entropies"])
close(fid)
number_of_bits = 16
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(log2.(guessing_entropies)), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="Logarithmic guessing entropy", ylabel="ChaCha state", xlabel="Operation number", clim=(0,15), dpi=300)
savefig(p, "./plots/heatmaps/32_16_bit_templates_LGE.png")
println("32: LGE: ", log2(mean(guessing_entropies)), " FSR: ", mean(success_rates))
