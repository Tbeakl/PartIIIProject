using HDF5, Plots
gr()

lower_bound_samples = 563000
upper_bound_samples = 587000

intermediate_value_index = 1
number_of_intermediate_values = 2800

intermediate_values_base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\intermediate_value_traces\\recording_profiling_"
traces_base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_profiling_"

trace_range_per_file = 0:249
file_range = 9:72

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
        sections_of_trace[(i-file_range[1])*length(trace_range_per_file)+j+1, :] = read(fid[string("power_", j)])[lower_bound_samples:(upper_bound_samples-1)]
    end
    close(fid)
end

# Need to get the correlations for each of the different parts of trace, will just try the hamming weight of the intermediate value for the time being

hamming_weight_values = Base.count_ones.(all_intermediate_values)
all_correlations = zeros(upper_bound_samples - lower_bound_samples, 4)
for i in 1:(upper_bound_samples-lower_bound_samples)
    for j in 1:4
        all_correlations[i, j] = cor(hamming_weight_values[:, j], sections_of_trace[:, i])^2
    end
end
p = plot(all_correlations, size=(1200, 800))
savefig(p, "./plots/hamming_weight_correlations.svg")

all_correlations_on_linear_model = zeros(upper_bound_samples - lower_bound_samples, 4)
for j in 1:4
    base_prediction_matrix = permutedims(hcat(digits.(all_intermediate_values[:, j], base=2, pad=9)...))
    base_prediction_matrix[:, end] .= 1
    for i in 1:(upper_bound_samples-lower_bound_samples)
        # Need to make the correct weights matrix for the prediction
        β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * sections_of_trace[:, i]
        predicted_values = base_prediction_matrix * β
        all_correlations_on_linear_model[i, j] = cor(predicted_values, sections_of_trace[:, i])^2
    end
end
p = plot(all_correlations_on_linear_model, size=(1200, 800))
savefig(p, "./plots/linear_models_correlations.svg")

# These correlation graphs seem to show that it is good to start each clock cycle at approximate 450/950 because that incleudes all the big peaks and is in a trough