using HDF5, Plots, Statistics, CUDA
plotly()

lower_bound_samples = 563000
upper_bound_samples = 587000

intermediate_value_index = 1
number_of_intermediate_values = 2800

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
sections_of_trace = zeros(Int16, length(trace_range_per_file) * length(file_range), upper_bound_samples - lower_bound_samples)

for i in file_range
    println(i)
    fid = h5open(string(intermediate_values_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        all_intermediate_values[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])
    end
    close(fid)

    fid = h5open(string(traces_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        power_trace = read(fid[string("power_", j)])
        difference_between_mean_and_power = argmin(power_trace) - mean_arg_min
        trimmed_power_trace = power_trace[50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
        sections_of_trace[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = trimmed_power_trace[lower_bound_samples:(upper_bound_samples-1)]
    end
    close(fid)
end
# These correlation graphs seem to show that it is good to start each clock cycle at approximate 450/950 because that incleudes all the big peaks and is in a trough
data_path = "D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_profiling.hdf5"
data_fid = h5open(data_path, "r")
intermediate_value_vector = read(data_fid["intermediate_values"])[1:16000, :]
dset = data_fid["downsampled_matrix"]
if HDF5.ismmappable(dset)
    dset = HDF5.readmmap(dset)
end
original_matrix_of_current_data = dset[1:16000, :]
close(data_fid)

sections_of_trace = CuArray(original_matrix_of_current_data)
all_intermediate_values = CuArray(intermediate_value_vector)

variances = var(sections_of_trace, dims=1)[1, :]

function calculate_NICV(traces, variances, all_intermediate_values, intermediate_index)
    intermediate_value_vector = all_intermediate_values[:, intermediate_index]
    all_mean_values = CUDA.zeros(size(traces)[2], 256)
    for j in 0:255
        all_mean_values[:, j + 1] = mean(traces[intermediate_value_vector .== j, :], dims=1)[1, :]
    end
    return var(all_mean_values, dims=2) ./ variances
end

all_NICV = CUDA.zeros(size(sections_of_trace)[2], 4)
for i in 1001:1004
    all_NICV[:, i - 1000] = calculate_NICV(sections_of_trace, variances, all_intermediate_values, i)
end
p = plot(Array(all_NICV))
p = plot(mean(Array(all_NICV), dims=2)[:, 1])
savefig(p, "./plots/NICVs.html")