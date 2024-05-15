# The first 1500 (files 1 to 6) will have there traces meaned to be compared to the other traces
# by checking that they are well correlated so that we can have confidence that the trigger has worked
using HDF5, StatsBase

samples_per_trace = 650_000

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

base_path = string(path_to_data, "captures/ChaChaRecordings_3/recording_profiling_")

summed_trace = zeros(Float32, samples_per_trace - 100)
trace_to_align_to = zeros(Int16, samples_per_trace)
all_gains::Vector{Float64} = []
all_offsets::Vector{Float64} = []


# Taken from https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal
@inline function allequal(x)
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i = 2:length(x)
        x[i] == e1 || return false
    end
    return true
end

# Should look to align all the traces correctly so the mean is correctly aligned with them
for i in 0:1
    fid = h5open(string(base_path, i, ".hdf5"))
    for j in 0:999
        println(i, " ", j)
        if i + j == 0
            trace_to_align_to = read(fid[string("power_", j)])
        else
            # Align the result to the mean trace
            raw_trace = read(fid[string("power_", j)])
            base_difference_between_mean_and_power = (argmin(raw_trace) - argmin(trace_to_align_to))
            lags_to_try = (-5:5) .+ base_difference_between_mean_and_power
            difference_between_mean_and_power = lags_to_try[argmax(crosscor(trace_to_align_to, raw_trace, lags_to_try))]
            trimmed_raw_trace = raw_trace[begin+50+difference_between_mean_and_power:end-(50-difference_between_mean_and_power)]
            summed_trace[:] .+= trimmed_raw_trace
        end
        push!(all_gains, read(fid[string("power_", j)]["gain"]))
        push!(all_offsets, read(fid[string("power_", j)]["offset"]))
    end
    close(fid)
end
summed_trace .+ trace_to_align_to[begin+50:end-50]
fid = h5open(string(path_to_data, "attack_profiling/32_volatile/mean_trace.hdf5"), "w")
@assert allequal(all_gains)
@assert allequal(all_offsets)
fid["mean_trace"] = summed_trace ./ 2000
attributes(fid["mean_trace"])["offset"] = all_offsets[1]
attributes(fid["mean_trace"])["gain"] = all_gains[1]
close(fid)