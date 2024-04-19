using HDF5, StatsBase, Statistics, Base.Threads

clock_cycle_sample_number = 405
number_of_intermediate_values = 2800
number_of_clock_cycles = (749500 ÷ 500)

intermediate_values_base_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/intermediate_value_traces/recording_profiling_"
traces_base_path = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/captures/ChaChaRecordings/recording_profiling_"

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/mean_trace.hdf5", "r")
mean_trace = read(fid["mean_trace"])
mean_arg_min = argmin(mean_trace)
close(fid)

trace_range_per_file = 0:249
# file_range = 73:328
file_range = 329:332
# file_range = 181:

number_of_samples_to_average_over = 20

all_intermediate_values = zeros(UInt8, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)

downsampled_matrix_20 = zeros(Float32, length(trace_range_per_file) * length(file_range), 1499)
# downsampled_matrix_100 = zeros(Float32, length(trace_range_per_file) * length(file_range), 7495)
# downsampled_matrix_250 = zeros(Float32, length(trace_range_per_file) * length(file_range), 2998)

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
        downsampled_trace_20 = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, 500)))
        # downsampled_trace_100 = collect(Iterators.map(mean, Iterators.partition(downsampled_trace_20, 5)))
        # downsampled_trace_250 = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, 250)))
        downsampled_matrix_20[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_20
        # downsampled_matrix_100[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_100
        # downsampled_matrix_250[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_250
    end
    close(fid)
end

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_500_traces_validation.hdf5", "w")
fid["intermediate_values"] = all_intermediate_values
fid["downsampled_matrix"] = downsampled_matrix_20
close(fid)

# fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_100_traces_validation.hdf5", "w")
# fid["intermediate_values"] = all_intermediate_values
# fid["downsampled_matrix"] = downsampled_matrix_100
# close(fid)

# fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/downsampled_250_traces_validation.hdf5", "w")
# fid["intermediate_values"] = all_intermediate_values
# fid["downsampled_matrix"] = downsampled_matrix_250
# close(fid)

# trace_range_per_file = 0:249
# file_range = 329:332

# all_intermediate_values = zeros(UInt8, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
# mean_value_per_cycle = zeros(Float64, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)
# min_value_per_cycle = zeros(Float32, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)
# max_value_per_cycle = zeros(Float32, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)

# number_of_samples_to_average_over = 10
# number_of_downsampled_samples_per_clock_cycle = 500 ÷ number_of_samples_to_average_over

# for i in file_range
#     println(i)
#     fid = h5open(string(intermediate_values_base_path, i, ".hdf5"))
#     for j in trace_range_per_file
#         all_intermediate_values[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])
#     end
#     close(fid)

#     fid = h5open(string(traces_base_path, i, ".hdf5"))
#     for j in trace_range_per_file
#         raw_trace = read(fid[string("power_", j)])
#         difference_between_mean_and_power = argmin(raw_trace) - mean_arg_min
#         trimmed_raw_trace = raw_trace[50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
#         trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(500-clock_cycle_sample_number)-1)]
#         downsampled_trace = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, number_of_samples_to_average_over)))
#         mean_value_per_cycle[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = collect(Iterators.map(mean, Iterators.partition(downsampled_trace, number_of_downsampled_samples_per_clock_cycle)))
#         min_value_per_cycle[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = collect(Iterators.map(minimum, Iterators.partition(downsampled_trace, number_of_downsampled_samples_per_clock_cycle)))
#         max_value_per_cycle[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = collect(Iterators.map(maximum, Iterators.partition(downsampled_trace, number_of_downsampled_samples_per_clock_cycle)))
#     end
#     close(fid)
# end

# data_fid = h5open("D:\\ChaChaData\\attack_profiling\\cycle_values_validation.hdf5", "w")
# data_fid["intermediate_values"] = all_intermediate_values
# data_fid["mean_value_per_cycle"] = mean_value_per_cycle
# data_fid["min_value_per_cycle"] = min_value_per_cycle
# data_fid["max_value_per_cycle"] = max_value_per_cycle
# close(data_fid)