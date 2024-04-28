using HDF5

path_to_results = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/constant_counter_4_8/"

initial_ranks = zeros(1000)
final_ranks = zeros(1000)

for i in 1:1000
    fid = h5open(string(path_to_results, i, ".hdf5"), "r")
    initial_ranks[i] = read(fid["initial_estimated_rank_log2"])
    final_ranks[i] = read(fid["final_estimated_rank_log2"])
    close(fid)
end

println("Average initial rank: ", mean(initial_ranks))
println("Min initial rank: ", minimum(initial_ranks))
println("Max initial rank: ", maximum(initial_ranks))
println("Average final rank: ", mean(final_ranks))
println("Min final rank: ", minimum(final_ranks))
println("Max final rank: ", maximum(final_ranks))