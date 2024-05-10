using HDF5, Base.Threads, StatsBase, Statistics, Random
include("../attacks/byte_template_attacks/template_attack_traces.jl")
include("../encryption/leakage_functions.jl")

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"
number_of_rounds_to_average_over = 100000
dimensions = 8

for signal_to_noise_ratio in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6]
    Random.seed!(1234)

    noise = noise_distribution_fixed_standard_dev(1.0, dimensions)
    mean_vectors = generate_mean_vectors(noise, signal_to_noise_ratio, 255)

    summed_guessing_entropy::Int64 = 0
    summed_correct_counts::Int64 = 0

    for _ in 1:number_of_rounds_to_average_over
        # Need to generate a random vector
        value = rand(0:255)
        new_vector = rand(noise) .+ mean_vectors[:, value + 1]
        permutation_of_results = sortperm(get_prob_dist_of_vector(mean_vectors', noise, new_vector); rev=true)
        cur_guess_entropy = findfirst(x -> x == value, (0:255)[permutation_of_results])
        summed_guessing_entropy += cur_guess_entropy
        summed_correct_counts += (cur_guess_entropy == 1)
    end
    println(signal_to_noise_ratio, " LGE: ", log2(summed_guessing_entropy / number_of_rounds_to_average_over), " SR: ", summed_correct_counts / number_of_rounds_to_average_over)
end