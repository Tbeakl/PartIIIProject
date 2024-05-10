using HDF5, Base.Threads, StatsBase, Statistics, Random, NPZ
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

template_set = "D"

mean_guessing_entropy::Float64 = 0
mean_success_rate::Float64 = 0

for i in 0:3
    # Need to read the file
    mean_guessing_entropy += mean(parse.(Float64, hcat(split.(replace.(readlines(
                path_to_data * "templates_Keccak/Result_Tables/Result_Tables/GE_table_" *
                template_set * lpad(i, 2, "0") * "_G0.txt")[13:2:end-2], "\\" => ""), "&")...)[2:end, :]))

    mean_success_rate += mean(parse.(Float64, hcat(split.(replace.(readlines(
                path_to_data * "templates_Keccak/Result_Tables/Result_Tables/SR_table_" *
                template_set * lpad(i, 2, "0") * "_G0.txt")[13:2:end-2], "\\" => ""), "&")...)[2:end, :]))
end

println(template_set,
    " SR: ", mean_success_rate / 4,
    " LGE: ", log2(mean_guessing_entropy / 4))