using HDF5, StatsBase,  Statistics

clock_cycle_sample_number = 405
number_of_intermediate_values = 2800
number_of_clock_cycles = (749500 ÷ 500)

intermediate_values_base_path = "D:\\ChaChaData\\intermediate_value_traces\\recording_profiling_"
traces_base_path = "D:\\ChaChaData\\captures\\ChaChaRecordings\\recording_profiling_"

fid = h5open("D:\\ChaChaData\\attack_profiling\\mean_trace.hdf5", "r")
mean_trace = read(fid["mean_trace"])
trimmed_mean_trace = mean_trace[50:(end-50)]
mean_arg_min = argmin(mean_trace)
close(fid)

trace_range_per_file = 0:249
file_range = 329:332

number_of_samples_to_average_over = 10
number_of_downsampled_samples_per_clock_cycle = 500 ÷ number_of_samples_to_average_over

all_intermediate_values = zeros(UInt8, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
downsampled_matrix = zeros(Float32, length(trace_range_per_file) * length(file_range), 74941)

for i in file_range
    println(i)
    fid = h5open(string(intermediate_values_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        all_intermediate_values[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])
    end
    close(fid)

    fid = h5open(string(traces_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        raw_trace = read(fid[string("power_", j)])
        difference_between_mean_and_power = argmin(raw_trace) - mean_arg_min
        trimmed_raw_trace = raw_trace[50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
        trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(500-clock_cycle_sample_number)-1)]
        downsampled_trace = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, number_of_samples_to_average_over)))
        downsampled_matrix[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace
    end
    close(fid)
end
fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_validation.hdf5", "w")
fid["intermediate_values"] = all_intermediate_values
fid["downsampled_matrix"] = downsampled_matrix
close(fid)
