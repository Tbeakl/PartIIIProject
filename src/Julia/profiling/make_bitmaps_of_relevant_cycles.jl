using HDF5, Plots, StatsBase, Statistics, Images

mean_trace_threshold_mean = 0.02
mean_trace_threshold_min = 0.02
mean_trace_threshold_max = 0.025

individual_trace_mean_NICV = 0.025
individual_trace_min_NICV = 0.025
individual_trace_max_NICV = 0.025

number_of_intermediate_values = 2800
number_of_clock_cycles = (749500 ÷ 500)

# Collect each 4 together because the operations are on 32 bits at a time
all_mean_NICV = zeros(number_of_clock_cycles, number_of_intermediate_values)
all_min_NICV = zeros(number_of_clock_cycles, number_of_intermediate_values)
all_max_NICV = zeros(number_of_clock_cycles, number_of_intermediate_values)

fid = h5open("D:\\ChaChaData\\attack_profiling\\NICV_aligned_traces.hdf5", "r")
for intermediate_value_index in 1:number_of_intermediate_values
    all_mean_NICV[:, intermediate_value_index] .= fid[string("mean_", intermediate_value_index)]
    all_min_NICV[:, intermediate_value_index] .= fid[string("min_", intermediate_value_index)]
    all_max_NICV[:, intermediate_value_index] .= fid[string("max_", intermediate_value_index)]
end
close(fid)

all_bitmasks = zeros(Bool, number_of_clock_cycles, number_of_intermediate_values)

for intermediate_value_index in 1:4:number_of_intermediate_values
    intermediate_values = intermediate_value_index .+ (0:3)
    mean_bitmask = StatsBase.mean(all_mean_NICV[:, intermediate_values], dims=2)[:, 1] .>= mean_trace_threshold_mean
    min_bitmask = StatsBase.mean(all_min_NICV[:, intermediate_values], dims=2)[:, 1] .>= mean_trace_threshold_min
    max_bitmask = StatsBase.mean(all_max_NICV[:, intermediate_values], dims=2)[:, 1] .>= mean_trace_threshold_max
    base_bit_mask = (mean_bitmask .| min_bitmask .| max_bitmask)
    for current_intermediate_index in intermediate_values
        mean_bitmask = all_mean_NICV[:, current_intermediate_index] .>= individual_trace_mean_NICV
        min_bitmask = all_min_NICV[:, current_intermediate_index] .>= individual_trace_min_NICV
        max_bitmask = all_max_NICV[:, current_intermediate_index] .>= individual_trace_max_NICV

        all_bitmasks[:, current_intermediate_index] = (base_bit_mask .| mean_bitmask .| min_bitmask .| max_bitmask)
    end
end

base_bitmasks = all_bitmasks

number_of_cycles_to_expand = 8

all_bitmasks = dilate(base_bitmasks; dims=1, r=number_of_cycles_to_expand)

fid = h5open("D:\\ChaChaData\\attack_profiling\\clock_cycles_bitmasks.hdf5", "w")
for intermediate_value_index in 1:number_of_intermediate_values
    fid[string("bitmask_", intermediate_value_index)] = all_bitmasks[:, intermediate_value_index]
end
close(fid)