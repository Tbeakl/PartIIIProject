# The first 1500 (files 1 to 6) will have there traces meaned to be compared to the other traces
# by checking that they are well correlated so that we can have confidence that the trigger has worked
using HDF5, StatsBase

samples_per_trace = 750_000
number_of_samples_to_average_over = 10
base_path = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_profiling_"

summed_trace = zeros(Float32, samples_per_trace ÷ number_of_samples_to_average_over)
all_gains::Vector{Float64} = []
all_offsets::Vector{Float64} = []


# Taken from https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal
@inline function allequal(x)
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end


for i in 1:6
    fid = h5open(string(base_path, i, ".hdf5"))
    for j in 0:249
        println(i, " ", j)
        summed_trace[:] .+= collect(Iterators.map(mean, Iterators.partition(read(fid[string("power_", j)]), number_of_samples_to_average_over)))
        push!(all_gains, read(fid[string("power_", j)]["gain"]))
        push!(all_offsets, read(fid[string("power_", j)]["offset"]))
    end
    close(fid)
end

all_gains

fid = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\mean_trace_downsampled.hdf5", "w")
@assert allequal(all_gains)
@assert allequal(all_offsets)

fid["mean_trace"] = summed_trace ./ 1500
attributes(fid["mean_trace"])["offset"] = all_offsets[1]
attributes(fid["mean_trace"])["gain"] = all_gains[1]
close(fid)