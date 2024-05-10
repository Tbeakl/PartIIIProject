using HDF5, StatsBase, Statistics, Base.Threads, Plots
plotly()
include("common_functions.jl")
clock_cycle_sample_number = 46
number_of_intermediate_values = 700
number_of_samples_per_clock_cycle = 50
number_of_clock_cycles = (759500 ÷ number_of_samples_per_clock_cycle)

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

intermediate_values_base_path = string(path_to_data, "intermediate_value_traces_8_on_32/recording_profiling_")
traces_base_path = string(path_to_data, "captures/ChaChaRecordings_8_on_32/recording_profiling_")

fid = h5open(string(path_to_data, "attack_profiling/8_on_32/mean_trace.hdf5"), "r")
mean_trace = read(fid["mean_trace"])
mean_arg_min = argmin(mean_trace)
close(fid)

trace_range_per_file = 0:999
# file_range = 46:57
file_range = 46:46
# file_range = 46:57
# file_range = 58:58
# file_range = 26:35

all_intermediate_values = zeros(UInt32, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)

# downsampled_matrix_20 = zeros(Float32, length(trace_range_per_file) * length(file_range), 1499)
# downsampled_matrix_50 = zeros(Float32, length(trace_range_per_file) * length(file_range), 759851)
downsampled_matrix_100 = zeros(Float32, length(trace_range_per_file) * length(file_range), 399925 * 2)
# downsampled_matrix_250 = zeros(Float32, length(trace_range_per_file) * length(file_range), 151971)
# downsampled_matrix_500 = zeros(Float32, length(trace_range_per_file) * length(file_range), 15198)

for i in file_range
    println(i)
    fid = h5open(string(intermediate_values_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        all_intermediate_values[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])
    end
    close(fid)

    fid = h5open(string(traces_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        println(j)
        trimmed_raw_trace = make_power_trace_trimmed_and_aligned_to_mean(mean_trace, read(fid[string("power_", j)]))
        trimmed_raw_trace = trimmed_raw_trace[clock_cycle_sample_number:(end-(number_of_samples_per_clock_cycle-clock_cycle_sample_number)-1)]
        # downsampled_trace_50 = trimmed_raw_trace #collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, 50)))
        downsampled_trace_100 = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, 1)))
        # downsampled_trace_250 = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, 5)))
        # downsampled_trace_500 = collect(Iterators.map(mean, Iterators.partition(trimmed_raw_trace, 50)))
        # downsampled_matrix_50[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_50
        downsampled_matrix_100[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_100
        # downsampled_matrix_250[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_250
        # downsampled_matrix_500[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = downsampled_trace_500
    end
    close(fid)
end

p = plot(size=(1500,1500))
plot!(p, downsampled_matrix_100[1, :])
plot!(p, downsampled_matrix_100[2, :])
plot!(p, downsampled_matrix_100[3, :])
plot!(p, downsampled_matrix_100[4, :])
plot!(p, downsampled_matrix_100[5, :])
plot!(p, trimmed_mean_trace)
# fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/8_on_32_trace_set/profiling_1_1.hdf5", "w")
# fid["intermediate_values"] = all_intermediate_values
# fid["downsampled_matrix"] = downsampled_matrix_50
# close(fid)

# fid = h5open(string(path_to_data, "attack_profiling/8_on_32/profiling_2_4.hdf5"), "w")
# fid["intermediate_values"] = all_intermediate_values
# fid["downsampled_matrix"] = downsampled_matrix_100
# close(fid)

# fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/8_on_32_trace_set/profiling_5_1.hdf5", "w")
# fid["intermediate_values"] = all_intermediate_values
# fid["downsampled_matrix"] = downsampled_matrix_250
# close(fid)

# fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/8_on_32_trace_set/profiling_50_1.hdf5", "w")
# fid["intermediate_values"] = all_intermediate_values
# fid["downsampled_matrix"] = downsampled_matrix_500
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