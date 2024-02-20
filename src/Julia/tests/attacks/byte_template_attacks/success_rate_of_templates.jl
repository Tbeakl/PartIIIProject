using NPZ, Distributions
include("../../../attacks/byte_template_attacks/template_attack_traces.jl")
include("../../../encryption/leakage_functions.jl")


function ASCON_key_success_rates()
    base_path_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_expect_b"
    base_path_noise = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_scov_b"
    # noise = noise_distribution_fixed_standard_dev(1., 8)

    expected_success_rates = [0.859 , 0.809 , 0.852 , 0.733 , 0.804 , 0.684 , 0.751 , 0.619 , 0.791 , 0.758 , 0.820 , 0.766 , 0.868 , 0.777 , 0.758 , 0.647]
    our_succes_rates = zeros(16)
    number_of_iterations = 10000

    for j in 0:15
        mean_vectors = transpose(npzread(string(base_path_mean_vectors, lpad(string(j % 16), 3, "0"), ".npy")))
        # Basically seems to make no difference to the values which are returned
        noise = noise_distribution_given_covaraince_matrix(npzread(string(base_path_noise, lpad(string(j % 16), 3, "0"), ".npy")))
        number_of_successes = 0
        for _ in 1:number_of_iterations
            # Generate a random value between 0 and 255 and put that in and see if we have guesses it correctly
            value_input = rand(0:255)
            prob_dist = make_prob_dist_for_byte(mean_vectors, noise, value_input)
            most_likely_value = findmax(prob_dist)[2] - 1
            if value_input == most_likely_value
                number_of_successes += 1
            end
        end
        our_succes_rates[j + 1] = number_of_successes / number_of_iterations
    end

    println(expected_success_rates)
    println(our_succes_rates)
end

function seeing_difference_expected_and_sampled_dist()
    value_to_sample = 0
    template_number = 0
    base_path_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_expect_b"
    noise = noise_distribution_fixed_standard_dev(1., 8)
    mean_vectors = transpose(npzread(string(base_path_mean_vectors, lpad(string(template_number % 16), 3, "0"), ".npy")))

    sampling_likelihoods = zeros(256)

    expected_likelihoods = pdf(noise, mean_vectors .- mean_vectors[:, value_to_sample + 1])


    number_of_iterations = 1000
    for _ in 1:number_of_iterations
        # Generate a random value between 0 and 255 and put that in and see if we have guesses it correctly
        likelihood_dist = make_prob_dist_for_byte(mean_vectors, noise, value_to_sample)
        sampling_likelihoods += likelihood_dist
    end

    sampling_likelihoods /= sum(sampling_likelihoods)
    expected_likelihoods /= sum(expected_likelihoods)

    println("Expected ", expected_likelihoods)
    println("Sampling ", sampling_likelihoods)

    println(sum(abs.(sampling_likelihoods .- expected_likelihoods)))
end

ASCON_key_success_rates()

attempt_at_unit_test()