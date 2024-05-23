using HDF5, Base.Threads, StatsBase, Statistics, Random, NPZ
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/evaluation/template_qualities/eight_bit_fragments_validation.hdf5"

fid = h5open(path_to_data, "r")
# mean_success_rate = mean(read(fid["intermediate_success_rate"]))
# mean_guessing_entropy = mean(read(fid["intermediate_guessing_entropies"]))
mean_success_rate = mean(read(fid["success_rates"]))
mean_guessing_entropy = mean(read(fid["guessing_entropies"]))
close(fid)

println("SR: ", mean_success_rate,
    " LGE: ", log2(mean_guessing_entropy))