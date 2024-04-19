using HDF5, StatsBase, Statistics, Base.Threads

function dilate_infront(bit_vector::AbstractVector, dilation_amount::Int64)
    if dilation_amount == 0
        return bit_vector
    end
    bit_vector = copy(bit_vector)
    end_locations_of_dilations = findall(x -> x == 1, diff(bit_vector))
    for loc in end_locations_of_dilations
        bit_vector[max(loc - dilation_amount + 1, 1):loc] .= true
    end
    return bit_vector
end

function dilate_after(bit_vector::AbstractVector, dilation_amount::Int64)
    if dilation_amount == 0
        return bit_vector
    end
    bit_vector = copy(bit_vector)
    start_locations_of_dilations = findall(x -> x == -1, diff(bit_vector))
    end_index = length(bit_vector)
    for loc in start_locations_of_dilations
        bit_vector[loc:min(loc + dilation_amount, end_index)] .= true
    end
    return bit_vector
end

clock_cycle_sample_number = 405
number_of_intermediate_values = 4
base_intermediate_value_index = 1
number_of_clock_cycles = (749500 ÷ 500)

intermediate_values_base_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/intermediate_value_traces/recording_profiling_"
traces_base_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/captures/ChaChaRecordings/recording_profiling_"
bitmask_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation.hdf5"


fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/mean_trace.hdf5", "r")
mean_trace = read(fid["mean_trace"])
mean_arg_min = argmin(mean_trace)
close(fid)

trace_range_per_file = 0:249
# file_range = 9:72
# file_range = 73:328
file_range = 329:332

number_of_samples_to_average_over = 25
number_of_downsampled_samples_per_clock_cycle = 500 ÷ number_of_samples_to_average_over

bitmask_fid = h5open(bitmask_path, "r")
cycle_bitmask = dilate_after(read(bitmask_fid[string("bitmask_", base_intermediate_value_index)]), 3)
cycle_bitmask = dilate_infront(cycle_bitmask, 4)
# cycle_bitmask = read(bitmask_fid[string("bitmask_", intermediate_value_index)])
close(bitmask_fid)

println(sum(cycle_bitmask))

sample_bitmask = repeat(cycle_bitmask, inner=500)[1:749401]

all_intermediate_values = zeros(UInt8, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
downsampled_matrix = zeros(Int16, length(trace_range_per_file) * length(file_range), sum(sample_bitmask))
for i in file_range
    println(i)
    fid = h5open(string(intermediate_values_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        all_intermediate_values[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])[base_intermediate_value_index:(base_intermediate_value_index + number_of_intermediate_values - 1)]
    end
    close(fid)

    fid = h5open(string(traces_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        raw_trace = read(fid[string("power_", j)])
        difference_between_mean_and_power = argmin(raw_trace) - mean_arg_min
        trimmed_raw_trace = raw_trace[50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
        trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(500-clock_cycle_sample_number)-1)]
        downsampled_matrix[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = trimmed_raw_trace[sample_bitmask]
    end
    close(fid)
end
fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/validation_1.hdf5", "w")
fid["intermediate_values"] = all_intermediate_values
fid["downsampled_matrix"] = downsampled_matrix
close(fid)