using Plots, HDF5, DSP, StatsBase, Statistics, DataStructures, GR
# plotly()
# gr()
usecolorscheme(2)
file_to_plot = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_profiling_156.hdf5"
fid = h5open(file_to_plot, "r")

fid_mean = h5open("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\attack_profiling\\mean_trace.hdf5", "r")
mean_trace = read(fid_mean["mean_trace"])
mean_offset = read(fid_mean["mean_trace"]["offset"])
mean_gain = read(fid_mean["mean_trace"]["gain"])
close(fid_mean)

trace_num_to_plot = 204

# power_trace = read(fid[string("trigger_", trace_num_to_plot)])
# power_offset = read(fid[string("trigger_", trace_num_to_plot)]["offset"])
# power_gain = read(fid[string("trigger_", trace_num_to_plot)]["gain"])

power_trace = read(fid[string("power_", trace_num_to_plot)])
power_offset = read(fid[string("power_", trace_num_to_plot)]["offset"])
power_gain = read(fid[string("power_", trace_num_to_plot)]["gain"])

# power_trace = collect(Iterators.map(mean, Iterators.partition(power_trace, number_of_samples_to_average_over)))
# power_trace = read(fid[string("power_", trace_num_to_plot, "_1")])
# power_offset = read(fid[string("power_", trace_num_to_plot, "_1")]["offset"])
# power_gain = read(fid[string("power_", trace_num_to_plot, "_1")]["gain"])
# trigger_trace = read(fid[string("trigger_", trace_num_to_plot)])

# difference_between_mean_and_power = argmin(power_trace) - argmin(mean_trace)
# trimmed_power_trace = power_trace[50 + difference_between_mean_and_power:end-(50 - difference_between_mean_and_power)]
# trimmed_mean_trace = mean_trace[50:(end-50)]

# p = plot([(trimmed_power_trace .* power_gain) .+ power_offset], label="Current trace", size=(1200,800), dpi=500)
# plot!(p, [(trimmed_mean_trace .* mean_gain) .+ mean_offset], label="Mean trace")

p = shade((power_trace .* power_gain) .+ power_offset, colormap=GR.COLORMAP_BLUESCALE, size=(1500, 500), dpi=300)
# savefig(p, "testing.png")
# savefig(p, "testing_plot.png")

# lags = -5:5
# best_lags::Vector{Int64} = []
# correlations::Vector{Float64} = []
# base_power_signal = read(fid["power_0"])
# for i in 1:249
#     power_sig = read(fid[string("power_", i)])
#     push!(best_lags, lags[argmax(crosscor(base_power_signal, power_sig, lags))])
#     push!(correlations, cor(base_power_signal, power_sig))
# end
# println(counter(best_lags))
# println(minimum(correlations))
close(fid)
p