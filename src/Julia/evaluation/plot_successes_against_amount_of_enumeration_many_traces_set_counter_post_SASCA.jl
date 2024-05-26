using Plots, HDF5, LaTeXStrings
include("../encryption/key_enumeration.jl")
gr()
base_path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

tickfontsize=11
guidefontsize=14
titlefontsize=16
legendfontsize=11

# Traces which have NaNs in them 347, 384, 403 need to potentially look in more detail at why they do and 
# if there is someway they can be improved
base_paths_to_counts::Vector{String} = base_path_to_data .* [
    "evaluation/set_counter_mean_10_8/",
    "evaluation/set_counter_prod_10_8/",
    "evaluation/set_counter_mean_10_16/",
    "evaluation/set_counter_prod_10_16/",
    "evaluation/random_counter_32_volatile_set_mean_10_8/",
    "evaluation/random_counter_32_volatile_set_prod_10_8/"
]

paths_to_actual_keys::Vector{String} = base_path_to_data .* "captures/" .* [
    "ChaChaRecordings_2/recording_attack_counter_constant_",
    "ChaChaRecordings_2/recording_attack_counter_constant_",
    "ChaChaRecordings_2/recording_attack_counter_constant_",
    "ChaChaRecordings_2/recording_attack_counter_constant_",
    "ChaChaRecordings_3/recording_attack_counter_constant_",
    "ChaChaRecordings_3/recording_attack_counter_constant_",
]

final_ranks::Vector{Vector{Number}} = []

labels::Vector{String} = [
"8-bit fragment mean",
"8-bit fragment product",
"16-bit fragment mean",
"16-bit fragment product",
"8-bit fragment volatile mean",
"8-bit fragment voaltile product"]

for (i, base_path_to_counts) in enumerate(base_paths_to_counts)
    current_final_ranks::Vector{Number} = []
    for trace_number in 1:1000
        if ispath(string(base_path_to_counts, trace_number, ".hdf5"))
            fid = h5open(string(base_path_to_counts, trace_number, ".hdf5"), "r")
            # Need to do all of the different types of real key enumeration, possibly need to increase
            # the number of iterations done on the unknown output version because it may not be correct
            current_estimated_rank = read(fid[string("final_estimated_rank_log2")])
            if current_estimated_rank < 20
                # Actually perform the key enumeration to find the solution
                # need to also get the key for this to know the correct values which is not stored in it 
                prob_tables = read(fid[string("final_likelihood_tables")])
                # Also need to read in the correct key so need to pick the correct file and trace within it
                # println(trace_number)
                file_number = (trace_number ÷ 100)
                trace_number_in_file = (trace_number - 1) % 100
                key_fid = h5open(paths_to_actual_keys[i] * string(file_number) * ".hdf5", "r")
                key = UInt32.(read(key_fid["power_"*string(trace_number_in_file)*"_0"]["key"]))
                close(key_fid)

                sorted_likelihood_matrix, sorted_values = make_matrices_of_values(eachrow(prob_tables))
                key_by_cluster = turn_key_into_cluster_values(key, 8)
                actual_rank = key_enumerate(sorted_likelihood_matrix, sorted_values, key_by_cluster, 1 << 20)
                if actual_rank >= (1 << 20)
                    println(base_path_to_counts, ": ", trace_number)
                    println(calculate_log_likelihood_of_key(key, eachrow(prob_tables), 8))
                    println(current_estimated_rank)
                end
                push!(current_final_ranks, log2(actual_rank))
            else
                push!(current_final_ranks, current_estimated_rank)
            end
            # if isnan(read(fid["entropy_over_time"])[end])
            #     println(base_path_to_counts, " ", trace_number)
            # end
            close(fid)
        end
    end
    push!(final_ranks, sort(current_final_ranks))
end

proportion = collect((1:1000) ./ 1000)

push!(proportion, 1.0)
for i in eachindex(base_paths_to_counts)
    push!(final_ranks[i], 256)
end

# Also add zero at the start of the trace and prepend the other intial rank to get the long line at the bottom
# and another zero zero at the very start
prepend!(proportion, [0, 0])
for i in eachindex(base_paths_to_counts)
    prepend!(final_ranks[i], [0, final_ranks[i][begin]])
end

# I think it would be good to have the products as dashed lines with the mean ones solid and have them as the same
# colour as each other
cur_colors = get_color_palette(:auto, plot_color(:white))

p = plot(size=(1500, 500),
    title="Proportion of keys successfully found after differing amounts of key enumeration",
    ylabel="Proportion",
    xlabel="Estimated number of keys required to be enumerated (log scale)",
    leftmargin=8Plots.mm,
    bottom_margin=8Plots.mm,
    legend=:outerright, legendcolumns=1, xlim=(0, 256), ylim=(0, 1),
    xticks=([0:32:256;], latexstring.("2^{" .* string.(0:32:256) .* "}")),
    yticks=([0:0.2:1;], latexstring.(0:0.2:1)),
    xtickfontsize=tickfontsize, 
    ytickfontsize=tickfontsize, 
    legendfontsize=legendfontsize,
    guidefontsize=guidefontsize)

for i in 1:2:length(base_paths_to_counts)
    cur_colors = get_color_palette(:auto, plot_color(:white))
    plot!(p, final_ranks[i], proportion, label=labels[i], linewidth=2, c=cur_colors[(i ÷ 2) + 1], linestyle=:solid)
    plot!(p, final_ranks[i + 1], proportion, label=labels[i + 1], linewidth=2, c=cur_colors[(i ÷ 2) + 1], linestyle=:dash)
end
p
savefig(p, "./plots/evaluation/real_attack_multi_trace_set_counter_post_SASCA.pdf")