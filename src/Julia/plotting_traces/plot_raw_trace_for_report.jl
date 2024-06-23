using Plots, HDF5, DSP, StatsBase, Statistics, DataStructures, GR, GRUtils # Plots, 
include("../profiling/common_functions.jl")
gr()
usecolorscheme(1)

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

file_to_plot = string(path_to_data, "captures/ChaChaRecordings_2/recording_attack_counter_constant_0.hdf5")
trace_num_to_plot = 1

fid = h5open(file_to_plot, "r")
power_trace = read(fid[string("power_", trace_num_to_plot, "_0")])
power_offset = read(fid[string("power_", trace_num_to_plot, "_0")]["offset"])
power_gain = read(fid[string("power_", trace_num_to_plot, "_0")]["gain"])
close(fid)

times = (0:749999) ./ 2_500_000_000

p = Plots.plot(times, power_trace .* power_gain .+ power_offset, size=(1500,500), ylabel="V", xlabel="Time (s)",
 xlim=(times[begin],times[end]), bottom_margin=10Plots.mm, left_margin=8Plots.mm, right_margin=6Plots.mm, label=false,
 guidefontsize=14, legendfontsize=14, titlefontsize=16, tickfontsize=10, dpi=300)
vspan!(p, [42_000, 532_000] ./ 2_500_000_000, label="Main rounds", alpha=.3)
vspan!(p, [550_000, 630_000] ./ 2_500_000_000, label="Add initial state", alpha=.3)
vspan!(p, [630_000, 730_000] ./ 2_500_000_000, label="Xor plaintext", alpha=.3)
Plots.savefig("./plots/raw_trace/raw_trace_different_parts_highlighted.svg")
Plots.savefig("./plots/raw_trace/raw_trace_different_parts_highlighted.pdf")
Plots.savefig("./plots/raw_trace/raw_trace_different_parts_highlighted.png")

# x = collect(1:length(power_trace)) ./ 2_500_000_00
y = (power_trace .* power_gain) .+ power_offset
p = GRUtils.Figure((1500,500))
p = GRUtils.shade!(p, times, y; xlabel="Time (s)", ylabel="V", colormap=-GR.COLORMAP_BLUESCALE, figsize=(1500,500))
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