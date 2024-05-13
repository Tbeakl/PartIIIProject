using HDF5, Plots, StatsBase, Statistics, Images
include("../chacha_factor_graph/heatmap_visualisation_intermediate_values.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")

individual_trace_mean_NICV = 0.004

number_of_intermediate_values = 700
number_of_bits = 8
templates_per_intermediate_value = 32 ÷ number_of_bits
number_of_clock_cycles = 15997 #(759500 ÷ 50)

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

# Collect each 4 together because the operations are on 32 bits at a time
all_mean_NICV = zeros(number_of_clock_cycles, number_of_intermediate_values * templates_per_intermediate_value)

for intermediate_value_index in 1:number_of_intermediate_values
    fid = h5open(string(path_to_data * "attack_profiling/8_on_32/correlations/", intermediate_value_index, ".hdf5"), "r")
    cur_inter_corrlations = read(fid["correlations"])
    for cluster_num in 1:templates_per_intermediate_value
        println(intermediate_value_index, " ", cluster_num)
        all_mean_NICV[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] = cur_inter_corrlations[cluster_num, :]
    end
    close(fid)
end

all_bitmasks = zeros(Bool, number_of_clock_cycles, number_of_intermediate_values * templates_per_intermediate_value)

for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        all_bitmasks[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] = all_mean_NICV[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] .>= individual_trace_mean_NICV
    end
end

mapping = turn_intermediate_name_to_intermediate_index(number_of_bits_per_template)
base_matrix = map(x -> x[begin:end-2], make_positions_to_var_names(number_of_bits_per_template, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
heatmap_of_counts = map_to_values.(mapped_matrix, Ref(sum(all_bitmasks, dims=1)[1, :]), Ref(32 ÷ number_of_bits_per_template))
p = heatmap(heatmap_of_counts, title="Number of interesting clock cycles", ylabel="Location in state by bytes", dpi=300)
savefig(p, "./plots/heatmaps/correlation_interesting_cycles_8_bits.png")


fid = h5open(path_to_data * "attack_profiling/8_on_32/clock_cycles_bitmasks.hdf5", "w")
for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        println(intermediate_value_index, " ", cluster_num)
        fid[string("bitmask_", intermediate_value_index, "_", cluster_num)] = all_bitmasks[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num]
    end
end
close(fid)