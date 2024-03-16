using Plots, HDF5, DSP, StatsBase, Statistics, DataStructures#, GR
plotly()
# gr()
# usecolorscheme(2)
file_to_plot = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings\\recording_attack_counter_from_one_1.hdf5"
fid = h5open(file_to_plot, "r")

trace_num_to_plot = 44

# power_trace = read(fid[string("trigger_", trace_num_to_plot)])
# power_offset = read(fid[string("trigger_", trace_num_to_plot)]["offset"])
# power_gain = read(fid[string("trigger_", trace_num_to_plot)]["gain"])

power_trace = read(fid[string("power_", trace_num_to_plot, "_5")])
power_offset = read(fid[string("power_", trace_num_to_plot, "_5")]["offset"])
power_gain = read(fid[string("power_", trace_num_to_plot, "_5")]["gain"])

# power_trace = read(fid[string("power_", trace_num_to_plot, "_1")])
# power_offset = read(fid[string("power_", trace_num_to_plot, "_1")]["offset"])
# power_gain = read(fid[string("power_", trace_num_to_plot, "_1")]["gain"])
# trigger_trace = read(fid[string("trigger_", trace_num_to_plot)])

p = plot([(power_trace .* power_gain) .+ power_offset], size=(1500,500), dpi=500)

# p = shade((power_trace .* power_gain) .+ power_offset, colormap=GR.COLORMAP_BLUESCALE, size=(1500, 500), dpi=300)
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