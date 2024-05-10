using HDF5, Base.Threads, StatsBase, Statistics, Random, NPZ
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"
number_of_rounds_to_average_over = 1000

template_set = "D"
number_of_templates::Int64 = 40

Random.seed!(1234)

summed_guessing_entropy::Int64 = 0
summed_correct_counts::Int64 = 0

for _ in 1:number_of_rounds_to_average_over
    # Need to pick a random template to use for the mean vectors
    template_set_number = lpad(rand(0:3), 2, "0")
    template_number = lpad(rand(0:(number_of_templates - 1)), 3, "0")

    mean_vectors = npzread(path_to_data *
                           "templates_Keccak/templateLDA_B_ID/template_" *
                           template_set * template_set_number *
                           "/template_expect_b" * template_number * ".npy")'

    noise = noise_distribution_given_covaraince_matrix(npzread(path_to_data *
                                                               "templates_Keccak/templateLDA_B_ID/template_" *
                                                               template_set * template_set_number *
                                                               "/template_scov_b" * template_number * ".npy"))

    # Need to generate a random vector
    value = rand(0:255)
    new_vector = rand(noise) .+ mean_vectors[:, value+1]
    permutation_of_results = sortperm(get_prob_dist_of_vector(mean_vectors', noise, new_vector); rev=true)
    cur_guess_entropy = findfirst(x -> x == value, (0:255)[permutation_of_results])
    summed_guessing_entropy += cur_guess_entropy
    summed_correct_counts += (cur_guess_entropy == 1)
end

println(template_set, " SR: ", summed_correct_counts / number_of_rounds_to_average_over, " LGE: ", log2(summed_guessing_entropy / number_of_rounds_to_average_over))