using NPZ, Distributions
include("../../../attacks/byte_template_attacks/template_attack_traces.jl")
include("../../../encryption/leakage_functions.jl")

function ASCON_key_success_rates()
    base_path_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_expect_b"
    base_path_noise = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_scov_b"
    # noise = noise_distribution_fixed_standard_dev(1., 8)

    expected_success_rates = [0.859 , 0.809 , 0.852 , 0.733 , 0.804 , 0.684 , 0.751 , 0.619 , 0.791 , 0.758 , 0.820 , 0.766 , 0.868 , 0.777 , 0.758 , 0.647]
    our_succes_rates = zeros(16)
    number_of_iterations = 1000

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

function Keccack_A_success_rates()
    base_path_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_A00\\template_expect_b"
    # base_path_noise = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_scov_b"
    noise = noise_distribution_fixed_standard_dev(1., 8)

    # expected_success_rates = [0.859 , 0.809 , 0.852 , 0.733 , 0.804 , 0.684 , 0.751 , 0.619 , 0.791 , 0.758 , 0.820 , 0.766 , 0.868 , 0.777 , 0.758 , 0.647]
    our_succes_rates = zeros(200)
    number_of_iterations = 1000

    for j in 0:199
        mean_vectors = transpose(npzread(string(base_path_mean_vectors, lpad(string(j % 200), 3, "0"), ".npy")))
        # Basically seems to make no difference to the values which are returned
        # noise = noise_distribution_given_covaraince_matrix(npzread(string(base_path_noise, lpad(string(j % 200), 3, "0"), ".npy")))
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

    # println(expected_success_rates)
    println(our_succes_rates)
end


ASCON_key_success_rates()

Keccack_A_success_rates()