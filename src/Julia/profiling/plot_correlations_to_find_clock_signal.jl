using HDF5, Plots, Statistics
plotly()

lower_bound_samples = 563000
upper_bound_samples = 587000

intermediate_value_index = 1
number_of_intermediate_values = 2800

intermediate_values_base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\intermediate_value_traces\\recording_profiling_"
traces_base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_profiling_"

fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\mean_trace.hdf5", "r")
mean_trace = read(fid["mean_trace"])
trimmed_mean_trace = mean_trace[50:(end-50)]
mean_arg_min = argmin(mean_trace)
close(fid)

trace_range_per_file = 0:249
file_range = 9:72 #9:72

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

# Need to get the correlations for each of the different parts of trace, will just try the hamming weight of the intermediate value for the time being

# hamming_weight_values = Base.count_ones.(all_intermediate_values)
# all_correlations = zeros(upper_bound_samples - lower_bound_samples, 4)
# for j in 1:4
#     weight_vector = hamming_weight_values[:, 4*j-1]
#     for i in 1:(upper_bound_samples-lower_bound_samples)
#         all_correlations[i, j] = cor(weight_vector, sections_of_trace[:, i])^2
#     end
# end
# p = plot(all_correlations, size=(1200, 800))
# savefig(p, "./plots/hamming_weight_correlations.html")

# all_correlations_on_linear_model = zeros(upper_bound_samples - lower_bound_samples, 4)
# for j in 1:4
#     base_prediction_matrix = permutedims(hcat(digits.(all_intermediate_values[:, j], base=2, pad=9)...))
#     base_prediction_matrix[:, end] .= 1
#     for i in 1:(upper_bound_samples-lower_bound_samples)
#         # Need to make the correct weights matrix for the prediction
#         β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * sections_of_trace[:, i]
#         predicted_values = base_prediction_matrix * β
#         all_correlations_on_linear_model[i, j] = cor(predicted_values, sections_of_trace[:, i])^2
#     end
# end
# p = plot(all_correlations_on_linear_model, size=(1200, 800))
# savefig(p, "./plots/linear_models_correlations.html")

# These correlation graphs seem to show that it is good to start each clock cycle at approximate 450/950 because that incleudes all the big peaks and is in a trough

value_columns = eachcol(sections_of_trace)
variances = var.(value_columns)

function calculate_normalised_interclass_variances(intermediate_value_vector, trace_values, variance)
    get_mean(j) = mean(trace_values[intermediate_value_vector .== j])
    return var(get_mean.(0:255)) / variance
end

all_NICV = zeros(upper_bound_samples - lower_bound_samples, 4)
for j in 1:4
    intermediate_value_vector = all_intermediate_values[:, j]
    all_NICV[:, j] = calculate_normalised_interclass_variances.(Ref(intermediate_value_vector), value_columns, variances)
end
p = plot(all_NICV, size=(1200, 800))
savefig(p, "./plots/NICVs.html")