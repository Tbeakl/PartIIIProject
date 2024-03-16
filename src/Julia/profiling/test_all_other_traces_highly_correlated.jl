using HDF5

base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_attack_counter_from_one_"
correlation_threshold = 0.98

fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\mean_trace.hdf5", "r")
mean_trace = read(fid["mean_trace"])
close(fid)

all_traces_to_exclude::Vector{String} = []
for file_number in 1:3
    fid = h5open(string(base_path, file_number, ".hdf5"), "r")
    # for i in 0:249
    #     corelation = cor(mean_trace, read(fid[string("power_", i)]))
    #     println(file_number, " ", i, ": ", corelation)
    #     if corelation <= correlation_threshold
    #         push!(all_traces_to_exclude, string(file_number, "_", i))
    #     end
    # end
    for i in 0:99
        for j in 0:9
            corelation = cor(mean_trace, read(fid[string("power_", i, "_", j)]))
            println(file_number, " ", i, " ", j, ": ", corelation)
            if corelation <= correlation_threshold
                push!(all_traces_to_exclude, string(file_number, "_", i, "_", j))
            end
        end
    end
    close(fid)
end

# fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\profiling_traces_less_than_0_98_correlation.hdf5", "w")
# fid["all_traces"] = all_traces_to_exclude
# close(fid)

# Done correlations on with triggers and attack traces and they are all above the threshold