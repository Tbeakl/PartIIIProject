using Plots, HDF5
gr()
base_path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

base_paths_to_counts::Vector{String} = base_path_to_data .* [
    "evaluation/attack_8_on_32_known/",
    "evaluation/attack_8_on_32_unknown/",
    "evaluation/random_counter_1_8/",
    "evaluation/random_counter_1_16/",
    "evaluation/random_counter_incremented_10_8/",
    "evaluation/random_counter_incremented_10_16/",
    "evaluation/incremented_counter_prod_10_8/",
    "evaluation/incremented_counter_prod_10_16/",
    "evaluation/set_counter_mean_10_8/",
    "evaluation/set_counter_mean_10_16/",
    "evaluation/set_counter_prod_10_8/",
    "evaluation/set_counter_prod_10_16/"]

final_ranks::Vector{Vector{Number}} = []

labels::Vector{String} = ["8-bit implementation\nsingle trace",
    "8-bit implementation\nunknown counter, nonce, output",
    "8-bit fragment\nsingle trace",
    "16-bit fragment\nsingle trace",
    "8-bit fragment 10 trace\nchanged counter mean",
    "16-bit fragment 10 trace\nchanged counter mean",
    "8-bit fragment 10 trace\nchanged counter product",
    "16-bit fragment 10 trace\nchanged counter product",
    "8-bit fragment 10 trace\nset counter mean",
    "16-bit fragment 10 trace\nset counter mean",
    "8-bit fragment 10 trace\nset counter product",
    "16-bit fragment 10 trace\nset counter product"]

for base_path_to_counts in base_paths_to_counts
    current_final_ranks::Vector{Number} = []
    for i in 1:1000
        if ispath(string(base_path_to_counts, i, ".hdf5"))
            fid = h5open(string(base_path_to_counts, i, ".hdf5"), "r")
            # Need to do all of the different types of real key enumeration, possibly need to increase
            # the number of iterations done on the unknown output version because it may not be correct
            current_estimated_rank = read(fid[string("final_estimated_rank_log2")])
            if current_estimated_rank < (1 << 20)
                # Actually perform the key enumeration to find the solution
                # need to also get the key for this to know the correct values which is not stored in it 
                # prob_tables = read(fid[string("")])
            end
            if isnan(read(fid["entropy_over_time"])[end])
                println(base_path_to_counts, " ", i)
            end
            push!(current_final_ranks, current_estimated_rank)
            close(fid)
        end
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
savefig(p, "./plots/evaluation/real_attacks_post_SASCA.pdf")
p