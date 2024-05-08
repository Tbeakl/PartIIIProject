using Plots, HDF5

base_paths_to_counts::Vector{String} = ["D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/8_on_32_random_counter_1_8/",
    "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/random_counter_1_8/",
    "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/random_counter_1_16/",]

final_ranks::Vector{Vector{Number}} = []

labels::Vector{String} = ["8-bit implementation", "8-bit fragments", "16-bit fragments marginalised\nto 8-bits"]

for base_path_to_counts in base_paths_to_counts
    current_final_ranks::Vector{Number} = []
    for i in 1:1000
        fid = h5open(string(base_path_to_counts, i, ".hdf5"), "r")
        push!(current_final_ranks, read(fid[string("initial_estimated_rank_log2")]))
        close(fid)
    end
    push!(final_ranks, sort(current_final_ranks))
end

p = plot(size=(1500, 500),
    title="Proportion of keys successfully found after differing amounts of key enumeration",
    ylabel="Proportion",
    xlabel="Log base 2 of estimated number of keys required to be enumerated",
    leftmargin=8Plots.mm,
    bottom_margin=6Plots.mm,
    legend=:outerright, legendcolumns=1, xlim=(0, 256), ylim=(0, 1))
proportion = (1:1000) ./ 1000

for i in eachindex(base_paths_to_counts)
    cur_colors = get_color_palette(:auto, plot_color(:white))
    plot!(p, final_ranks[i], proportion, label=labels[i])
end
savefig(p, "./plots/evaluation/real_single_trace_pre_SASCA.pdf")
p