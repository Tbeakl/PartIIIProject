using Plots, HDF5

base_path_to_counts = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/8_on_32_random_counter_unknown_output_counter_nonce_1_8/"

initial_ranks::Vector{Number} = []
final_ranks::Vector{Number} = []

for i in 1:1000
    if i != 280
        fid = h5open(string(base_path_to_counts, i, ".hdf5"), "r")
        push!(initial_ranks, read(fid[string("initial_estimated_rank_log2")]))
        push!(final_ranks, read(fid[string("final_estimated_rank_log2")]))
        close(fid)
    end
end

sort!(initial_ranks)
sort!(final_ranks)

p = plot(size=(1500, 500),
    title="Proportion of keys successfully found after differing amounts of key enumeration",
    ylabel="Proportion",
    xlabel="Log base 2 of estimated number of keys required to be enumerated",
    leftmargin=8Plots.mm,
    bottom_margin=6Plots.mm,
    legend=:outerright, legendcolumns=1, xlim=(0, 256), ylim=(0, 1))
proportion = (1:999) ./ 999
cur_colors = get_color_palette(:auto, plot_color(:white))
plot!(p, initial_ranks, proportion, label=string("Pre-SASCA"))
plot!(p, final_ranks, proportion, label=string("Post-SASCA"))
savefig(p, "./plots/evaluation/real_8_on_32_unknown_output.pdf")
p