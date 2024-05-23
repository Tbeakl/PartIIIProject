using Plots, HDF5

# Make plots of varying speeds of heatmap to reflect the time which is taken by each of the different methods

forwards_backwards_fps::Int64 = 5
simple_fps::Int64 = 10
ends_fps::Int64 = 50

x_label_name::String = "Operation number"
y_label_name::String = "ChaCha state"
path_to_data::String = "C:/Users/henry/Documents/PartIIIProject/data/"


# Load in the different visualisations which can be done
function make_heatmap(additional_path::String, fps::Int64, number_of_bits::Int64, gif_path::String)
    visualisation_of_entropies::Vector{Matrix{Float64}} = []
    fid = h5open(path_to_data * additional_path, "r")
    i::Int64 = 0
    while haskey(fid, string("entropies_", i))
        push!(visualisation_of_entropies, read(fid[string("entropies_", i)]))
        i += 1
    end
    close(fid)
    anim = @animate for i in eachindex(visualisation_of_entropies)
        heatmap(visualisation_of_entropies[i];
            title=string("Entropies of variables after ", i - 1, " iterations"),
            clim=(0, number_of_bits),
            xlabel=x_label_name,
            ylabel=y_label_name,
            yticks=([0.5:(32 ÷ number_of_bits):((512 ÷ number_of_bits) + 1);], 0:16))
    end
    p = gif(anim, gif_path, fps=fps)
    return p
end

p = make_heatmap("evaluation/heatmap_data/schedule_changes/real_8bit_2_clusters_known_forwards_backwards.hdf5",
    forwards_backwards_fps,
    2,
    "./plots/schedule_animations/real_8bit_2_clusters_known_forwards_backwards.gif")

p = make_heatmap("evaluation/heatmap_data/schedule_changes/real_8bit_2_clusters_known_all.hdf5",
    simple_fps,
    2,
    "./plots/schedule_animations/real_8bit_2_clusters_known_all.gif")

p = make_heatmap("evaluation/heatmap_data/schedule_changes/real_8bit_2_clusters_known_ends.hdf5",
    ends_fps,
    2,
    "./plots/schedule_animations/real_8bit_2_clusters_known_ends.gif")

p = make_heatmap("evaluation/heatmap_data/schedule_changes/real_8bit_2_clusters_forwards_backwards.hdf5",
    forwards_backwards_fps,
    2,
    "./plots/schedule_animations/real_8bit_2_clusters_forwards_backwards.gif")

p = make_heatmap("evaluation/heatmap_data/schedule_changes/real_8bit_2_clusters_all.hdf5",
    simple_fps,
    2,
    "./plots/schedule_animations/real_8bit_2_clusters_all.gif")

p = make_heatmap("evaluation/heatmap_data/schedule_changes/real_8bit_2_clusters_ends.hdf5",
    ends_fps,
    2,
    "./plots/schedule_animations/real_8bit_2_clusters_ends.gif")