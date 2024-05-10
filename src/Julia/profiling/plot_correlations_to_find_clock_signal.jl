using HDF5, Plots, Statistics, LaTeXStrings
include("../encryption/leakage_functions.jl")
include("common_functions.jl")
# gr()
plotly()

lower_bound_samples = 750_000 #563000
upper_bound_samples = 799_900 #587000

intermediate_value_index = 700
number_of_intermediate_values = 700
path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

intermediate_values_base_path = string(path_to_data, "intermediate_value_traces_8_on_32/recording_profiling_")
traces_base_path = string(path_to_data, "captures/ChaChaRecordings_8_on_32/recording_profiling_")

fid = h5open(string(path_to_data, "attack_profiling/mean_trace_8_on_32.hdf5"), "r")
mean_trace = read(fid["mean_trace"])
trimmed_mean_trace = mean_trace[50:(end-50)]
mean_arg_min = argmin(mean_trace)
close(fid)

trace_range_per_file = 0:999
file_range = 2:9

all_intermediate_values = zeros(UInt32, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
sections_of_trace = zeros(Int16, length(trace_range_per_file) * length(file_range), upper_bound_samples - lower_bound_samples)

Threads.@threads for i in file_range
    println(i)
    fid = h5open(string(intermediate_values_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        all_intermediate_values[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])
    end
    close(fid)

    fid = h5open(string(traces_base_path, i, ".hdf5"))
    for j in trace_range_per_file
        trimmed_power_trace = make_power_trace_trimmed_and_aligned_to_mean(mean_trace, read(fid[string("power_", j)]))
        sections_of_trace[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = trimmed_power_trace[lower_bound_samples:(upper_bound_samples-1)]
    end
    close(fid)
end
# These correlation graphs seem to show that it is good to start each clock cycle at approximate 450/950 because that incleudes all the big peaks and is in a trough
# data_path = "D:\\ChaChaData\\attack_profiling\\downsampled_10_traces_profiling.hdf5"
# data_fid = h5open(data_path, "r")
# intermediate_value_vector = read(data_fid["intermediate_values"])[1:16000, :]
# dset = data_fid["downsampled_matrix"]
# if HDF5.ismmappable(dset)
#     dset = HDF5.readmmap(dset)
# end
# original_matrix_of_current_data = dset[1:16000, :]
# close(data_fid)

# sections_of_trace = original_matrix_of_current_data # CuArray(original_matrix_of_current_data)
# all_intermediate_values = intermediate_value_vector # CuArray(intermediate_value_vector)

byte_intermediate_values = hcat(byte_values_for_input.(all_intermediate_values[:, intermediate_value_index])...)'

variances = var(sections_of_trace, dims=1)[1, :]

function calculate_NICV(traces, variances, byte_intermediate_values, intermediate_index)
    intermediate_value_vector = byte_intermediate_values[:, intermediate_index]
    all_mean_values = zeros(size(traces)[2], 256)
    for j in 0:255
        all_mean_values[:, j+1] = mean(traces[intermediate_value_vector.==j, :], dims=1)[1, :]
    end
    return var(all_mean_values, dims=2) ./ variances
end

all_NICV = zeros(size(sections_of_trace)[2], 4)
for i in 1:4
    all_NICV[:, i] = calculate_NICV(sections_of_trace, variances, byte_intermediate_values, i)
end

hamming_weight_values = Base.count_ones.(byte_intermediate_values)
all_correlations = zeros(upper_bound_samples - lower_bound_samples, 4)
for j in 1:4
    weight_vector = hamming_weight_values[:, j]
    Threads.@threads for i in 1:(upper_bound_samples-lower_bound_samples)
        all_correlations[i, j] = cor(weight_vector, sections_of_trace[:, i])^2
    end
end

all_correlations_on_linear_model = zeros(upper_bound_samples - lower_bound_samples, 4)
for j in 1:4
    base_prediction_matrix = permutedims(hcat(digits.(byte_intermediate_values[:, j], base=2, pad=9)...))
    base_prediction_matrix[:, end] .= 1
    Threads.@threads for i in 1:(upper_bound_samples-lower_bound_samples)
        # Need to make the correct weights matrix for the prediction
        β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * sections_of_trace[:, i]
        predicted_values = base_prediction_matrix * β
        all_correlations_on_linear_model[i, j] = cor(predicted_values, sections_of_trace[:, i])^2
    end
end

p = plot(all_NICV[:, 1], label=L"K_1", size=(1000, 300), ylabel="NICV", xlabel="Sample number", left_margin=5Plots.mm, bottom_margin=6Plots.mm, title="NICV for individual samples")
plot!(p, all_NICV[:, 2], label=L"K_2")
plot!(p, all_NICV[:, 3], label=L"K_3")
plot!(p, all_NICV[:, 4], label=L"K_4")
plot!(p, mean(all_NICV, dims=2)[:, 1], label=L"K_\mu")
# savefig(p, "./plots/clock_signal_detection/NICV_clock_signal_8_on_32.svg")
# savefig(p, "./plots/clock_signal_detection/NICV_clock_signal_8_on_32.pdf")
# savefig(p, "./plots/NICVs.html")

p = plot(all_correlations[:, 1], label=L"K_1", size=(1000, 300), ylabel=L"R^2", xlabel="Sample number", left_margin=5Plots.mm, bottom_margin=6Plots.mm, title="Hamming weight correlation for individual samples")
plot!(p, all_correlations[:, 2], label=L"K_2")
plot!(p, all_correlations[:, 3], label=L"K_3")
plot!(p, all_correlations[:, 4], label=L"K_4")
plot!(p, mean(all_correlations, dims=2)[:, 1], label=L"K_\mu")
# savefig(p, "./plots/clock_signal_detection/Hamming_weight_clock_signal_8_on_32.svg")
# savefig(p, "./plots/clock_signal_detection/Hamming_weight_clock_signal_8_on_32.pdf")

p = plot(all_correlations_on_linear_model[:, 1], label=L"K_1", size=(1000, 300), ylabel=L"R^2", xlabel="Sample number", left_margin=5Plots.mm, bottom_margin=6Plots.mm, title="Linear bit model correlation for individual samples")
plot!(p, all_correlations_on_linear_model[:, 2], label=L"K_2")
plot!(p, all_correlations_on_linear_model[:, 3], label=L"K_3")
plot!(p, all_correlations_on_linear_model[:, 4], label=L"K_4")
plot!(p, mean(all_correlations_on_linear_model, dims=2)[:, 1], label=L"K_\mu")
# savefig(p, "./plots/clock_signal_detection/Linear_model_clock_signal_8_on_32.svg")
# savefig(p, "./plots/clock_signal_detection/Linear_model_clock_signal_8_on_32.pdf")
