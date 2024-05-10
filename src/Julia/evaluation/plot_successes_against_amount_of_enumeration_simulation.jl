using Plots, HDF5


path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"
base_path_to_counts = path_to_data * "evaluation/simulation_unknown_output_nonce_counter/"

signal_to_noise_ratios = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]

initial_ranks::Dict{Number,Vector{Number}} = Dict()
final_ranks::Dict{Number,Vector{Number}} = Dict()

for snr in signal_to_noise_ratios
    current_initial_ranks::Vector{Float64} = []
    current_final_ranks::Vector{Float64} = []
    for i in 1:100
        fid = h5open(string(base_path_to_counts, snr, "/", i, ".hdf5"), "r")
        push!(current_initial_ranks, read(fid["initial_estimated_rank_log2"]))
        push!(current_final_ranks, read(fid["final_estimated_rank_log2"]))
        close(fid)
    end
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