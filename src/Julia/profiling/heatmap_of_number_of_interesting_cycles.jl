using HDF5, Plots
include("common_functions.jl")
include("../chacha_factor_graph/heatmap_visualisation_intermediate_values.jl")
include("../chacha_factor_graph/heatmap_visualisation.jl")

number_of_intermediate_values = 700
number_of_bits = 8
templates_per_intermediate_value = 32 ÷ number_of_bits
number_of_clock_cycles = 15198

all_bitmasks = zeros(UInt8, number_of_clock_cycles, number_of_intermediate_values * templates_per_intermediate_value)

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation_8_on_32.hdf5", "r")
for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        println(intermediate_value_index, " ", cluster_num)
        all_bitmasks[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] = read(fid[string("bitmask_", intermediate_value_index, "_", cluster_num)])
    end
end
close(fid)

counts_per_template = Int64.(sum(all_bitmasks, dims=1)[1,:])

mapping = turn_intermediate_name_to_intermediate_index(number_of_bits)
base_matrix = map(x->x[begin:end-2], make_positions_to_var_names(number_of_bits, 1)[1])
mapped_matrix = map(x -> mapping[x], base_matrix)
reduced_to_single_loc_per_intermediate_value = mapped_matrix[begin:end, :]
counts_of_interesting_clockcycles = map_to_values.(reduced_to_single_loc_per_intermediate_value, Ref(counts_per_template), Ref(4))

p = heatmap(counts_of_interesting_clockcycles, title="Number of interesting clock cycles", ylabel="Location in state", dpi=300)
savefig(p, "./plots/heatmaps/correlation_interesting_cycles_8_bits.png")