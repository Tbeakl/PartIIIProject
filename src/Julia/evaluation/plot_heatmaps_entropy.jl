using Plots, HDF5

path_to_data_for_heatmaps = "./data/evaluation/heatmap_data/"

# Read in the correct rounds for each in turn
fid = h5open(path_to_data_for_heatmaps * "loopy_adds/2.hdf5", "r")
visualisation_data = read(fid["entropies_50"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 50 iterations",
    clim=(0, 2),
    yticks=([0.5:16:257;], 0:16))
savefig(p, "./plots/factor_graph_findings/loopy_add_fifty_rounds.png")

fid = h5open(path_to_data_for_heatmaps * "tree_adds/2.hdf5", "r")
visualisation_data = read(fid["entropies_50"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 50 iterations",
    clim=(0, 2),
    yticks=([0.5:16:257;], 0:16))
savefig(p, "./plots/factor_graph_findings/tree_add_fifty_rounds.png")

fid = h5open(path_to_data_for_heatmaps * "tree_adds/1.hdf5", "r")
visualisation_data = read(fid["entropies_147"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 147 iterations",
    clim=(0, 1),
    yticks=([0.5:32:513;], 0:16))
savefig(p, "./plots/factor_graph_findings/chacha_1bit_before_end.png")

fid = h5open(path_to_data_for_heatmaps * "tree_adds/2.hdf5", "r")
visualisation_data = read(fid["entropies_91"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 91 iterations",
    clim=(0, 2),
    yticks=([0.5:16:257;], 0:16))
savefig(p, "./plots/factor_graph_findings/chacha_2bit_before_end.png")

fid = h5open(path_to_data_for_heatmaps * "tree_adds/4.hdf5", "r")
visualisation_data = read(fid["entropies_68"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 68 iterations",
    clim=(0, 4),
    yticks=([0.5:8:129;], 0:16))
savefig(p, "./plots/factor_graph_findings/chacha_4bit_before_end.png")

fid = h5open(path_to_data_for_heatmaps * "tree_adds/8.hdf5", "r")
visualisation_data = read(fid["entropies_61"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 61 iterations",
    clim=(0, 8),
    yticks=([0.5:4:65;], 0:16))
savefig(p, "./plots/factor_graph_findings/chacha_8bit_before_end.png")

fid = h5open(path_to_data_for_heatmaps * "xor_adds/2.hdf5", "r")
visualisation_data = read(fid["entropies_42"])
close(fid)
p = heatmap(visualisation_data, dpi=300, ylabel="ChaCha state",
    xlabel="Operation number",
    title="Entropies of variables after 42 iterations",
    clim=(0, 2),
    yticks=([0.5:16:257;], 0:16))
savefig(p, "./plots/factor_graph_findings/xor_replaced_add_round_42.png")