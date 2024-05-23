using Plots, HDF5
include("../chacha_factor_graph/heatmap_visualisation_intermediate_values.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")


all_bitmasks = zeros(UInt8, 2800, 15997)
fid = h5open("./data/attack_profiling/8_on_32/clock_cycles_bitmasks.hdf5", "r")
for intermediate_value_index in 1:700
    for cluster_num in 1:4
        all_bitmasks[4*(intermediate_value_index-1)+cluster_num, :] = read(fid[string("bitmask_", intermediate_value_index, "_", cluster_num)])
    end
end
close(fid)
number_of_bits = 8
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(sum(all_bitmasks, dims=2)[:, 1]), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="Number of interesting clock cycles",
 ylabel="ChaCha state", xlabel="Operation number",
  dpi=300,
  yticks=([0.5:(32 ÷ number_of_bits):(512 ÷ number_of_bits) + 1;], 0:16))
savefig(p, "./plots/heatmaps/correlation_interesting_cycles_8_bits.png")

all_bitmasks = zeros(UInt8, 2800, 6498)
fid = h5open("./data/attack_profiling/32_volatile/clock_cycles_bitmasks.hdf5", "r")
for intermediate_value_index in 1:700
    for cluster_num in 1:4
        all_bitmasks[4*(intermediate_value_index-1)+cluster_num, :] = read(fid[string("bitmask_", intermediate_value_index, "_", cluster_num)])
    end
end
close(fid)
number_of_bits = 32
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(sum(all_bitmasks, dims=2)[begin:4:end, 1]), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="Number of interesting clock cycles", ylabel="ChaCha state", xlabel="Operation number", dpi=300,
yticks=([0.5:(32 ÷ number_of_bits):(512 ÷ number_of_bits) + 1;], 0:16))
savefig(p, "./plots/heatmaps/correlation_interesting_cycles_32_volatile_bits.png")


all_bitmasks = zeros(UInt8, 2800, 1499)
fid = h5open("./data/attack_profiling/32/clock_cycles_bitmasks.hdf5", "r")
for intermediate_value_index in 1:700
    for cluster_num in 1:4
        all_bitmasks[4*(intermediate_value_index-1)+cluster_num, :] = read(fid[string("bitmask_", intermediate_value_index, "_", cluster_num)])
    end
end
close(fid)
number_of_bits = 32
mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(sum(all_bitmasks, dims=2)[begin:4:end, 1]), Ref(32 ÷ number_of_bits))
p = heatmap(heatmap_of_counts, title="Number of interesting clock cycles", ylabel="ChaCha state", xlabel="Operation number", dpi=300,
yticks=([0.5:(32 ÷ number_of_bits):(512 ÷ number_of_bits) + 1;], 0:16))
savefig(p, "./plots/heatmaps/correlation_interesting_cycles_32.png")