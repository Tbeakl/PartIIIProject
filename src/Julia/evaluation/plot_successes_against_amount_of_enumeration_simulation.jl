using Plots, HDF5, LaTeXStrings

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"
base_path_to_counts = path_to_data * "evaluation/simulation_known_output_nonce_counter/"

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
    xlabel="Estimated number of keys required to be enumerated (log scale)",
    leftmargin=8Plots.mm,
    bottom_margin=8Plots.mm,
    legend=:outerright, legendcolumns=1, xlim=(0, 256), ylim=(0, 1),
    xticks=([0:32:256;], latexstring.("2^{" .* string.(0:32:256) .* "}")),
    yticks=([0:0.2:1;], latexstring.(0:0.2:1)),
    xtickfont=font(10), 
    ytickfont=font(10), 
    legendfont=font(10))
proportion = collect((1:100) ./ 100)
cur_colors = get_color_palette(:auto, plot_color(:white))

# Add on a final value at the top left of each so that the lines continue along the top
push!(proportion, 1.0)
for i in eachindex(signal_to_noise_ratios)
    push!(initial_ranks[signal_to_noise_ratios[i]], 256)
    push!(final_ranks[signal_to_noise_ratios[i]], 256)
end
prepend!(proportion, [0, 0])
for i in eachindex(signal_to_noise_ratios)
    prepend!(initial_ranks[signal_to_noise_ratios[i]], [0, initial_ranks[signal_to_noise_ratios[i]][begin]])
    prepend!(final_ranks[signal_to_noise_ratios[i]], [0, final_ranks[signal_to_noise_ratios[i]][begin]])
end

for i in eachindex(signal_to_noise_ratios)
    plot!(p, initial_ranks[signal_to_noise_ratios[i]], proportion, c=cur_colors[i], linestyle=:dash, label=string("Pre-SASCA ", signal_to_noise_ratios[i], " SNR"), linewidth=2)
    plot!(p, final_ranks[signal_to_noise_ratios[i]], proportion, c=cur_colors[i], linestyle=:solid, label=string("Post-SASCA ", signal_to_noise_ratios[i], " SNR"), linewidth=2)
end
p
savefig(p, "./plots/evaluation/simulated_signal_to_noise_ratios_known_output.pdf")


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
    xlabel="Estimated number of keys required to be enumerated (log scale)",
    leftmargin=8Plots.mm,
    bottom_margin=8Plots.mm,
    legend=:outerright, legendcolumns=1, xlim=(0, 256), ylim=(0, 1),
    xticks=([0:32:256;], latexstring.("2^{" .* string.(0:32:256) .* "}")),
    yticks=([0:0.2:1;], latexstring.(0:0.2:1)),
    xtickfont=font(10), 
    ytickfont=font(10), 
    legendfont=font(10))
cur_colors = get_color_palette(:auto, plot_color(:white))
proportion = collect((1:100) ./ 100)
# Add on a final value at the top left of each so that the lines continue along the top
push!(proportion, 1.0)
for i in eachindex(signal_to_noise_ratios)
    push!(initial_ranks[signal_to_noise_ratios[i]], 256)
    push!(final_ranks[signal_to_noise_ratios[i]], 256)
end
prepend!(proportion, [0, 0])
for i in eachindex(signal_to_noise_ratios)
    prepend!(initial_ranks[signal_to_noise_ratios[i]], [0, initial_ranks[signal_to_noise_ratios[i]][begin]])
    prepend!(final_ranks[signal_to_noise_ratios[i]], [0, final_ranks[signal_to_noise_ratios[i]][begin]])
end

for i in eachindex(signal_to_noise_ratios)
    plot!(p, initial_ranks[signal_to_noise_ratios[i]], proportion, c=cur_colors[i], linestyle=:dash, label=string("Pre-SASCA ", signal_to_noise_ratios[i], " SNR"), linewidth=2)
    plot!(p, final_ranks[signal_to_noise_ratios[i]], proportion, c=cur_colors[i], linestyle=:solid, label=string("Post-SASCA ", signal_to_noise_ratios[i], " SNR"), linewidth=2)
end
p
savefig(p, "./plots/evaluation/simulated_signal_to_noise_ratios_unknown_output.pdf")