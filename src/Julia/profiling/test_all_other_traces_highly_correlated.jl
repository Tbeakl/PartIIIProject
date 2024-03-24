using HDF5

base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_attack_counter_from_one_"
correlation_threshold = 0.98
number_of_samples_to_average_over = 10

fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\mean_trace.hdf5", "r")
mean_trace = read(fid["mean_trace"])
trimmed_mean_trace = mean_trace[50:(end - 50)]
mean_arg_min = argmin(mean_trace)
close(fid)

all_traces_to_exclude::Vector{String} = []
all_correlations_excluded::Vector{Float64} = []
for file_number in 7:360
    fid = h5open(string(base_path, file_number, ".hdf5"), "r")
    # for i in 0:249
    #     # power_trace = collect(Iterators.map(mean, Iterators.partition(read(fid[string("power_", i)]), number_of_samples_to_average_over)))
    #     power_trace = read(fid[string("power_", i)])
    #     difference_between_mean_and_power = argmin(power_trace) - mean_arg_min
    #     trimmed_power_trace = power_trace[50 + difference_between_mean_and_power:end-(50 - difference_between_mean_and_power)]

    #     correlation = cor(trimmed_mean_trace, trimmed_power_trace)
    #     println(file_number, " ", i, ": ", correlation)
    #     if correlation <= correlation_threshold
    #         push!(all_traces_to_exclude, string(file_number, "_", i))
    #         push!(all_correlations_excluded, correlation)
    #     end
    # end
    for i in 0:99
        for j in 0:9
            power_trace = read(fid[string("power_", i, "_", j)])
            difference_between_mean_and_power = argmin(power_trace) - mean_arg_min
            trimmed_power_trace = power_trace[50 + difference_between_mean_and_power:end-(50 - difference_between_mean_and_power)]
            correlation = cor(trimmed_mean_trace, trimmed_power_trace)
            println(file_number, " ", i, " ", j, ": ", correlation)
            if correlation <= correlation_threshold
                push!(all_traces_to_exclude, string(file_number, "_", i, "_", j))
                push!(all_correlations_excluded, correlation)
            end
        end
    end
    close(fid)
end

# fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\profiling_traces_less_than_0_98_correlation_aligned_on_first_min.hdf5", "w")
# fid["all_traces"] = all_traces_to_exclude
# fid["all_correlations"] = all_correlations_excluded
# close(fid)

# Done correlations on with triggers and attack traces and they are all above the threshold

# With aligning on the minimum sample there is no trace which has a correlation below the threshold pretty much
# all are 0.998/997 some might be fracitonally lower