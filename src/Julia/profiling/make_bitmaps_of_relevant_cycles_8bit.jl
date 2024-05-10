using HDF5, Plots, StatsBase, Statistics, Images

individual_trace_mean_NICV = 0.004

number_of_intermediate_values = 700
number_of_bits = 8
templates_per_intermediate_value = 32 ÷ number_of_bits
number_of_clock_cycles = 15997 #(759500 ÷ 50)

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

# Collect each 4 together because the operations are on 32 bits at a time
all_mean_NICV = zeros(number_of_clock_cycles, number_of_intermediate_values * templates_per_intermediate_value)

fid = h5open(path_to_data * "attack_profiling/8_on_32/COR_aligned_traces.hdf5", "r")
for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        println(intermediate_value_index, " ", cluster_num)
        all_mean_NICV[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] = read(fid[string("mean_", intermediate_value_index, "_", cluster_num)])
    end
end
close(fid)

all_bitmasks = zeros(Bool, number_of_clock_cycles, number_of_intermediate_values * templates_per_intermediate_value)

for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        all_bitmasks[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] = all_mean_NICV[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] .>= individual_trace_mean_NICV
    end
end

# base_bitmasks = all_bitmasks

# interesting_clock_cycle_counts = sum(base_bitmasks, dims=1)[begin:4:end]

# key_counts = interesting_clock_cycle_counts[1:8]

# intermediate_value_cycles = interesting_clock_cycle_counts[29:end-16]

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
#     append!(intermediate_value_locations, [string(i, "_add") for i in 1:16])
#     return intermediate_value_locations
# end

# intermediate_locations = make_break_down_of_values()

# # # Need to break down the cycles into different ways of having counts
# open("counts_interesting_clock_cycles_8_on_32.csv", "w") do file
#     # First do the key
#     write(file, "Key\n")
#     for i in 5:12
#         write(file, string(i, ", ", key_counts[i-4], "\n"))
#     end
#     # First do the adds
#     write(file, "Add\n")
#     for i in 1:16
#         write(file, string(i, ", ", string(intermediate_value_cycles[intermediate_locations.==string(i, "_add")])[2:end-1], "\n"))
#     end
#     write(file, "Rotate\n")
#     for i in 1:16
#         if sum(intermediate_locations .== string(i, "_rot")) > 0
#             write(file, string(i, ", ", string(intermediate_value_cycles[intermediate_locations.==string(i, "_rot")])[2:end-1], "\n"))
#         end
#     end
# end

fid = h5open(path_to_data * "attack_profiling/8_on_32/clock_cycles_bitmasks.hdf5", "w")
for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        println(intermediate_value_index, " ", cluster_num)
        fid[string("bitmask_", intermediate_value_index, "_", cluster_num)] = all_bitmasks[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num]
    end
end
close(fid)