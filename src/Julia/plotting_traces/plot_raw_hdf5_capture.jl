using Plots, HDF5, DSP, StatsBase, Statistics, DataStructures, GR, GRUtils # Plots, 
include("../profiling/common_functions.jl")

plotly()
# gr()
usecolorscheme(1)

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

file_to_plot = string(path_to_data, "captures/ChaChaRecordings_8_on_32/recording_profiling_2.hdf5")
fid = h5open(file_to_plot, "r")

fid_mean = h5open(string(path_to_data, "attack_profiling/mean_trace_8_on_32.hdf5"), "r")
mean_trace = read(fid_mean["mean_trace"])
mean_offset = read(fid_mean["mean_trace"]["offset"])
mean_gain = read(fid_mean["mean_trace"]["gain"])
close(fid_mean)

trace_num_to_plot = 250

# power_trace = read(fid[string("trigger_", trace_num_to_plot)])
# power_offset = read(fid[string("trigger_", trace_num_to_plot)]["offset"])
# power_gain = read(fid[string("trigger_", trace_num_to_plot)]["gain"])

power_trace = read(fid[string("power_", trace_num_to_plot)])
power_offset = read(fid[string("power_", trace_num_to_plot)]["offset"])
power_gain = read(fid[string("power_", trace_num_to_plot)]["gain"])
close(fid)

# p = plot(power_trace .* power_gain .+ power_offset, size=(1500,500), ylabel="V", xlabel="Sample number",
#  xlim=(0,750_000), bottom_margin=10Plots.mm, left_margin=8Plots.mm, right_margin=6Plots.mm, label=false,
#  title="Raw trace")
# vspan!(p, [42_000, 532_000], label="Main rounds", alpha=.3)
# vspan!(p, [550_000, 630_000], label="Add initial state", alpha=.3)
# vspan!(p, [630_000, 730_000], label="Xor plaintext", alpha=.3)
# savefig("./plots/raw_trace_different_parts_highlighted.svg")
# savefig("./plots/raw_trace_different_parts_highlighted.pdf")
# power_trace = collect(Iterators.map(mean, Iterators.partition(power_trace, number_of_samples_to_average_over)))
# power_trace = read(fid[string("power_", trace_num_to_plot, "_1")])
# power_offset = read(fid[string("power_", trace_num_to_plot, "_1")]["offset"])
# power_gain = read(fid[string("power_", trace_num_to_plot, "_1")]["gain"])
# trigger_trace = read(fid[string("trigger_", trace_num_to_plot)])

trimmed_power_trace = make_power_trace_trimmed_and_aligned_to_mean(mean_trace, power_trace)

p = Plots.plot([(trimmed_power_trace .* power_gain) .+ power_offset], label="Current trace", size=(1200, 800), dpi=500)
Plots.plot!(p, [(mean_trace .* mean_gain) .+ mean_offset], label="Mean trace")
# p = plot(size=(1500,200))
# colorscheme(1)

# x = collect(1:length(power_trace)) ./ 2_500_000_00
# y = (power_trace .* power_gain) .+ power_offset
# p = GRUtils.Figure((1500,500))
# p = GRUtils.shade!(p, x, y; xlabel="Time (s)", ylabel="V", title="Raw trace", colormap=-GR.COLORMAP_BLUESCALE, figsize=(1500,500))
# GRUtils.savefig("test.png", p)
# p = shade(x, y, colormap=-GR.COLORMAP_BLUESCALE, size=(1500, 500), dpi=300, ylabel="V", xlabel="Time (s)",
# xlim=(0, 750_000), title="Raw trace")
# GRUtils.savefig("test.png", p)
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
p