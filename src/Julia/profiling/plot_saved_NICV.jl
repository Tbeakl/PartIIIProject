using Plots, HDF5, StatsBase, Statistics
plotly()
fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/NICV_aligned_traces_2.hdf5", "r")

number_of_intermediate_values = 700 #2800
number_of_clock_cycles = (749500 ÷ 500)

all_mean_NICV = zeros(number_of_clock_cycles, 4 * number_of_intermediate_values)

for intermediate_value_index in 1:number_of_intermediate_values
    for tem_num in 1:4
    all_mean_NICV[:, 4 * (intermediate_value_index - 1) + 3] .= fid[string("mean_", intermediate_value_index, "_", tem_num)]
    end
end
close(fid)

base_intermediate_value = 1

p = plot(size=(1200, 800))
plot!(p, all_mean_NICV[:, base_intermediate_value], label="Byte 1")
plot!(p, all_mean_NICV[:, base_intermediate_value + 1], label="Byte 2")
plot!(p, all_mean_NICV[:, base_intermediate_value + 2], label="Byte 3")
plot!(p, all_mean_NICV[:, base_intermediate_value + 3], label="Byte 4")
plot!(p, mean(all_mean_NICV[:, base_intermediate_value .+ (0:3)], dims=2)[:, 1], label="Mean")
hline!(p, [0.02])
hline!(p, [0.025])
# plot!(p, 0.2 .* base_bitmasks[:, base_intermediate_value][:, 1], label="Included points")
# plot!(p, 0.2 .* cycle_bitmask, label="Included points again")
