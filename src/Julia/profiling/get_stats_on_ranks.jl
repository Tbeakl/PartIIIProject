using HDF5

# path_to_results = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/8_on_32_random_counter_1_8/"
path_to_results = "D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/evaluation/8_on_32_random_counter_unknown_output_counter_nonce_1_8/"
initial_ranks::Vector{Real} = []
final_ranks::Vector{Real} = []

for i in 1:1000
    if ispath(string(path_to_results, i, ".hdf5"))
        fid = h5open(string(path_to_results, i, ".hdf5"), "r")
        push!(initial_ranks, read(fid["initial_estimated_rank_log2"]))
        push!(final_ranks, read(fid["final_estimated_rank_log2"]))
        close(fid)
    end
end

println("Average initial rank: ", mean(initial_ranks))
println("Min initial rank: ", minimum(initial_ranks))
println("Max initial rank: ", maximum(initial_ranks))
println("Average final rank: ", mean(final_ranks))
println("Min final rank: ", minimum(final_ranks))
println("Max final rank: ", maximum(final_ranks))