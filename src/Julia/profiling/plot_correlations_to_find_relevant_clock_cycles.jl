using HDF5, Plots, Base.Threads, StatsBase, CUDA, Statistics
plotly()
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
file_range = 9:9#9:72

all_intermediate_values = zeros(UInt8, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
mean_value_per_cycle = zeros(Float64, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)
min_value_per_cycle = zeros(Float32, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)
max_value_per_cycle = zeros(Float32, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)

number_of_samples_to_average_over = 10
number_of_downsampled_samples_per_clock_cycle = 500 ÷ number_of_samples_to_average_over

number_of_bits = 1
number_of_clusters = 8 ÷ number_of_bits

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
        mean_value_per_cycle[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = collect(Iterators.map(mean, Iterators.partition(downsampled_trace, number_of_downsampled_samples_per_clock_cycle)))
        min_value_per_cycle[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = collect(Iterators.map(minimum, Iterators.partition(downsampled_trace, number_of_downsampled_samples_per_clock_cycle)))
        max_value_per_cycle[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = collect(Iterators.map(maximum, Iterators.partition(downsampled_trace, number_of_downsampled_samples_per_clock_cycle)))
    end
    close(fid)
end

function calculate_NICV(traces, variances, all_intermediate_values, intermediate_index, number_of_bits, cluster_number)
    intermediate_value_vector = (all_intermediate_values[:, intermediate_index] .>> (number_of_bits * (cluster_number - 1))) .& ((1 << number_of_bits) - 1)
    all_mean_values = CUDA.zeros(size(traces)[2], 1 << number_of_bits)
    for j in 0:((1 << number_of_bits) - 1)
        all_mean_values[:, j + 1] = mean(traces[intermediate_value_vector .== j, :], dims=1)[1, :]
    end
    return var(all_mean_values, dims=2) ./ variances
end

all_intermediate_values = CuArray(all_intermediate_values)
mean_value_per_cycle = CuArray{Float32}(mean_value_per_cycle)
min_value_per_cycle = CuArray{Float32}(min_value_per_cycle)
max_value_per_cycle = CuArray{Float32}(max_value_per_cycle)

mean_variances = var(mean_value_per_cycle, dims=1)[1, :]
min_variances = var(min_value_per_cycle, dims=1)[1, :]
max_variances = var(max_value_per_cycle, dims=1)[1, :]

all_NICV_mean_value = CUDA.zeros(number_of_intermediate_values * number_of_clusters, size(max_value_per_cycle)[2])
all_NICV_max_value = CUDA.zeros(number_of_intermediate_values * number_of_clusters, size(max_value_per_cycle)[2])
all_NICV_min_value = CUDA.zeros(number_of_intermediate_values * number_of_clusters, length(max_value_per_cycle[1, :]))

for intermediate_value_index in 1005:1009
    for cluster_num in 1:number_of_clusters
        println(intermediate_value_index, " ", cluster_num)
        all_NICV_mean_value[(number_of_clusters * (intermediate_value_index - 1)) + cluster_num, :] = calculate_NICV(mean_value_per_cycle, mean_variances, all_intermediate_values, intermediate_value_index, number_of_bits, cluster_num)
        all_NICV_min_value[(number_of_clusters * (intermediate_value_index - 1)) + cluster_num, :] = calculate_NICV(min_value_per_cycle, min_variances, all_intermediate_values, intermediate_value_index, number_of_bits, cluster_num)
        all_NICV_max_value[(number_of_clusters * (intermediate_value_index - 1)) + cluster_num, :] = calculate_NICV(max_value_per_cycle, max_variances, all_intermediate_values, intermediate_value_index, number_of_bits, cluster_num)
    end
end

p = plot(mean(Array(all_NICV_mean_value)[8033:8064, :], dims=1)[1, :], legend=false)
p = plot(Array(all_NICV_mean_value)[8033:8064, :]', legend=false)
fid = h5open("D:\\ChaChaData\\attack_profiling\\clock_cycles_bitmasks.hdf5", "r")
cycle_bitmasks = zeros(4, size(max_value_per_cycle)[2])
for i in 1:4
    cycle_bitmasks[i, :] = read(fid[string("bitmask_", i + 1004)])
end
close(fid)
plot!(p, cycle_bitmasks' ./ 10, color=:red)

# fid = h5open("D:\\ChaChaData\\attack_profiling\\NICV_aligned_traces.hdf5", "w")
# for intermediate_value_index in 1:number_of_intermediate_values
#     println(intermediate_value_index)
#     fid[string("mean_", intermediate_value_index)] = Array(all_NICV_mean_value[intermediate_value_index, :])
#     fid[string("min_", intermediate_value_index)] = Array(all_NICV_min_value[intermediate_value_index, :])
#     fid[string("max_", intermediate_value_index)] = Array(all_NICV_max_value[intermediate_value_index, :])
# end
# close(fid)

# fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_traces_for_clock_cycle.hdf5", "w")
# fid["mean_values"] = Array(mean_value_per_cycle)
# fid["min_values"] = Array(min_value_per_cycle)
# fid["max_values"] = Array(max_value_per_cycle)
# fid["intermediate_values"] = Array(all_intermediate_values)
# close(fid)