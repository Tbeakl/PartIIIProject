using HDF5, Plots, Base.Threads, StatsBase
println(Threads.nthreads())
function main()
    clock_cycle_sample_number = 405
    number_of_intermediate_values = 2800
    number_of_clock_cycles = (749500 ÷ 500)

    intermediate_values_base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\intermediate_value_traces\\recording_profiling_"
    traces_base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_profiling_"

    fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\mean_trace.hdf5", "r")
    mean_trace = read(fid["mean_trace"])
    trimmed_mean_trace = mean_trace[50:(end-50)]
    mean_arg_min = argmin(mean_trace)
    close(fid)

    trace_range_per_file = 0:249
    file_range = 9:9 #9:72

    all_intermediate_values = zeros(UInt8, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
    mean_value_per_cycle = zeros(Float64, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)
    min_value_per_cycle = zeros(Float32, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)
    max_value_per_cycle = zeros(Float32, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)

    number_of_samples_to_average_over = 10
    number_of_downsampled_samples_per_clock_cycle = 500 ÷ number_of_samples_to_average_over

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

    # Need to work out which clock cycles relate to this intermediate value which requries defining a cut off for what is considered

    # Need to find the all the values which have already been created because we don't want to recreate them
    # fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/linear_model_correlations.hdf5", "cw")

    # close(fid)

    Threads.@threads for intermediate_value_index in 1:2800
        all_correlations_on_linear_model_mean_value = zeros(length(max_value_per_cycle[1, :]))
        all_correlations_on_linear_model_min_value = zeros(length(max_value_per_cycle[1, :]))
        all_correlations_on_linear_model_max_value = zeros(length(max_value_per_cycle[1, :]))
        println(intermediate_value_index)
        base_prediction_matrix = permutedims(hcat(digits.(all_intermediate_values[:, intermediate_value_index], base=2, pad=9)...))
        base_prediction_matrix[:, end] .= 1
        for i in 1:length(max_value_per_cycle[1, :])
            # Need to make the correct weights matrix for the prediction
            β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * mean_value_per_cycle[:, i]
            predicted_values = base_prediction_matrix * β
            all_correlations_on_linear_model_mean_value[i] = cor(predicted_values, mean_value_per_cycle[:, i])

            β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * min_value_per_cycle[:, i]
            predicted_values = base_prediction_matrix * β
            all_correlations_on_linear_model_min_value[i] = cor(predicted_values, min_value_per_cycle[:, i])

            β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * max_value_per_cycle[:, i]
            predicted_values = base_prediction_matrix * β
            all_correlations_on_linear_model_max_value[i] = cor(predicted_values, max_value_per_cycle[:, i])
        end
        fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/linear_model_correlations_aligned_traces.hdf5", "cw")
        fid[string("mean_", intermediate_value_index)] = all_correlations_on_linear_model_mean_value
        fid[string("min_", intermediate_value_index)] = all_correlations_on_linear_model_min_value
        fid[string("max_", intermediate_value_index)] = all_correlations_on_linear_model_max_value
        close(fid)
    end
    p = plot(all_correlations_on_linear_model_min_value, size=(1200, 800))
end

main()