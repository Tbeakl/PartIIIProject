using Plots, NPZ, Statistics, StatsBase

number_of_points_in_bin = 2
initial_trace = npzread("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\test_captures\\1_raw_capture.npy")

downsampled_power_trace = collect(Iterators.map(mean, Iterators.partition(initial_trace[2, :], number_of_points_in_bin)))
downsampled_power_trace ./= maximum(downsampled_power_trace)

# initial_trace = npzread("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\test_captures\\0_raw_capture.npy")
downsampled_trigger_trace = collect(Iterators.map(mean, Iterators.partition(initial_trace[1, :], number_of_points_in_bin)))
downsampled_trigger_trace ./= maximum(downsampled_trigger_trace)
# Just try reducing the number of poitns significantly by averaging together like 50 in a go

start_point = 125000
end_point = 175000
p = plot([downsampled_power_trace[start_point:end_point] downsampled_trigger_trace[start_point:end_point]], labels=["Power trace" "Trigger trace"], size=(1200,700))
savefig(p, "initial_plot.svg")

lags = -20:20
best_lags::Vector{Int64} = []


trigger_sig_1 = npzread("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\test_captures\\1_raw_capture.npy")[1, :]
for i in 1:49
    trigger_sig_2 = npzread(string("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\test_captures\\",i, "_raw_capture.npy"))[1, :]
    push!(best_lags, lags[argmax(crosscor(trigger_sig_1, trigger_sig_2, lags))])
end
p = plot()
for i in 0:49
    trigger_sig = npzread(string("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\test_captures\\",i, "_raw_capture.npy"))[1, :]
    plot!(p, trigger_sig[8000:9000])
end
plot!(p, size=(1000, 1000))

# By inspecting the trigger signal it appears that 1 is a very consistent point so it is probably a good
# threshold to use for the trigger

# Also need to look at cutting down the length of the recording because there is a significant amount after the end
# which is not actually useful because it is not doing any useful work

# Think it cas also be good to limit the recording lengths to 740,000 samples because the encryption has finished by that
# point so that means we can reduce the caputer length to like 1480 (maybe just 1500) cycles of recording.

# This means assuming 1500 * 500 * 2 = 1.5MB of storage per captured trace (plus a bit extra for the values for offset and gain etc.)
# Need to know how many traces will likely be required for storage but I think this should manage to fit into 
# my unix home of like 100GB.