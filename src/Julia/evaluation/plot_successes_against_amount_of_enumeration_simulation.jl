using Plots, HDF5

base_path_to_counts = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/simulation_unknown_output_nonce_counter/"

signal_to_noise_ratios = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]

initial_ranks::Dict{Number,Vector{Number}} = Dict()
final_ranks::Dict{Number,Vector{Number}} = Dict()

for snr in signal_to_noise_ratios
    current_initial_ranks::Vector{Float64} = []
    current_final_ranks::Vector{Float64} = []
    fid = h5open(string(base_path_to_counts, snr, "_successes.hdf5"), "r")
    for i in 1:100
        push!(current_initial_ranks, read(fid[string("initial_estimated_rank_log2_", i)]))
        push!(current_final_ranks, read(fid[string("final_estimated_rank_log2_", i)]))
    end
    close(fid)
    initial_ranks[snr] = sort(current_initial_ranks)
    final_ranks[snr] = sort(current_final_ranks)
end

p = plot(size=(1500, 500),
    title="Proportion of keys successfully found after differing amounts of key enumeration",
    ylabel="Proportion",
    xlabel="Log base 2 of estimated number of keys required to be enumerated",
    leftmargin=8Plots.mm,
    bottom_margin=6Plots.mm,
    legend=:outerright, legendcolumns=1, xlim=(0,256), ylim=(0,1))
proportion = (1:100) ./ 100
cur_colors = get_color_palette(:auto, plot_color(:white))
for i in eachindex(signal_to_noise_ratios)
    plot!(p, initial_ranks[signal_to_noise_ratios[i]], proportion, c=cur_colors[i], linestyle=:dash, label=string("Pre-SASCA ", signal_to_noise_ratios[i]))
    plot!(p, final_ranks[signal_to_noise_ratios[i]], proportion, c=cur_colors[i], linestyle=:solid, label=string("Post-SASCA ", signal_to_noise_ratios[i]))
end
savefig(p, "./plots/evaluation/simulated_signal_to_noise_ratios_unknown_output.pdf")