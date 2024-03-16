using Plots, Base.Threads, Random, NPZ, FileIO, Dates
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("template_attack_traces.jl")

number_of_bits = 2
dimensions = 8
initial_number_of_iterations = 250
number_of_iterations_of_ends = 0
rounds_for_ends = 5
number_of_encryption_traces = 1

noise = noise_distribution_fixed_standard_dev(1.0, dimensions)
base_path_key_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_A00\\template_expect_b"
base_path_intermediate_add_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_C00\\template_expect_b"
base_path_intermediate_rot_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_Keccak\\templateLDA_B_ID\\template_B00\\template_expect_b"

file_to_write_to_name = string("key_A_add_C_rot_B_runs_bits_2_", number_of_encryption_traces, ".csv")
fileFileIOStream = open(file_to_write_to_name, "a")

variables = Dict{String,Variable{Factor}}()
factors = Dict{String,Factor{Variable}}()
variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]

for encryption_run_number in 1:number_of_encryption_traces
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
    add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
end

additional_variables::Set{String} = Set{String}()
additional_factors::Set{String} = Set{String}()

if number_of_encryption_traces > 1
    add_equality_between_keys_and_nonces(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors)
    add_adds_between_counters(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors, additional_variables)
end

all_variables = []
all_factors = []

internal_factors = [union(factors_by_round[:]..., additional_factors...)...]
internal_variables = [union(variables_by_round[:]..., additional_variables...)...]

variables_at_ends = [union(variables_by_round[begin:rounds_for_ends]..., variables_by_round[end-rounds_for_ends-1:end]...)...]
factors_at_ends = [union(factors_by_round[begin:rounds_for_ends]..., factors_by_round[end-rounds_for_ends-1:end]...)...]

# Potentially could look at just the very start and end for checking for if we have reached convergence because we 
# don't really need to solve the centre of the graph if the start and end are correct and will not really change any more
number_of_end_rounds_to_check = 2
vars_to_check = [union(variables_by_round[begin:number_of_end_rounds_to_check]..., variables_by_round[end-number_of_end_rounds_to_check-1:end]...)...]
calculate_entropy_of_ends() = sum([variables[var].current_entropy for var in vars_to_check])

# This appears to be having some weird thing happening where the intial distributions 
# are not being added to the results because it seems to be completely independent amount of entropy
# reagrdless of the templates which have been chosen

for i in 1:25
    add_byte_key_template_to_variable = byte_template_path_to_function(base_path_key_mean_vectors, noise, 200)
    add_byte_intermediate_add_template_to_variable = byte_template_path_to_function(base_path_intermediate_add_mean_vectors, noise, 40)
    add_byte_intermediate_rot_template_to_variable = byte_template_path_to_function(base_path_intermediate_rot_mean_vectors, noise, 40)

    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
    adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]

    key = generate_random_key()
    nonce = generate_random_nonce()
    counter = generate_random_counter()

    for encryption_run_number in 1:number_of_encryption_traces
        encryption_trace = encrypt_collect_trace(key, nonce, counter, byte_values_for_input)
        encryption_output = encrypt(key, nonce, counter)
        location_execution_counts = zeros(Int64, 16)
        chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, encryption_run_number)
        add_starting_constant_values(variables, factors, number_of_bits, encryption_run_number)
        add_initial_nonce_and_counter_dist(variables, factors, number_of_bits, byte_values_for_input.(nonce), byte_values_for_input(counter), encryption_run_number, add_byte_key_template_to_variable)
        add_initial_key_dist(variables, factors, number_of_bits, byte_values_for_input.(key), encryption_run_number, add_byte_key_template_to_variable)
        # Need to add some noisy distributions to the initial values for the counters to see how that does

        for i in 1:16
            add_byte_key_template_to_variable(byte_values_for_input(encryption_output[i]), variables, factors, number_of_bits, string(i, "_", location_execution_counts[i]), encryption_run_number)
        end
        println("Starting to add the factors for the trace")
        add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, encryption_run_number, add_byte_intermediate_add_template_to_variable, add_byte_intermediate_rot_template_to_variable)
        println("Added the factors for the trace")
        counter = counter + 0x1
    end

    if number_of_encryption_traces > 1
        # Add equality constraint between all the encryption runs and the add constraint between the counters
        add_equality_between_keys_and_nonces(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors)
        add_adds_between_counters(variables, factors, number_of_bits, number_of_encryption_traces, additional_factors, additional_variables)
    end

    if length(all_factors) <= 0
        # Need to set all factors and all variables with the intial distributions included as ones to propagate about
        all_variables = [keys(variables)...]
        all_factors = [keys(factors)...]
    end

    # Do an initial push out of all the starting distributions 
    Threads.@threads for fact_name in all_factors
        factor_to_variable_messages(factors[fact_name])
    end

    Threads.@threads for var_name in all_variables
        variable_to_factor_messages(variables[var_name])
    end

    for add_nums in adds_by_round
        for add_num in add_nums
            for encryption_run in 1:number_of_encryption_traces
                belief_propagate_through_add(variables, factors, number_of_bits, add_num, encryption_run)
            end
        end
    end

    tot_entropy_over_time::Vector{Float64} = []
    update_all_entropies(variables, all_variables)
    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])

    for i in 1:initial_number_of_iterations
        println(i)
        Threads.@threads for var_name in internal_variables
            variable_to_factor_messages(variables[var_name], .8)
        end
        Threads.@threads for fact_name in internal_factors
            factor_to_variable_messages(factors[fact_name], .8)
        end
        update_all_entropies(variables, all_variables)
        push!(tot_entropy_over_time, total_entropy_of_graph(variables))
        println(tot_entropy_over_time[end])
        if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6
            break
        end
    end

    if number_of_iterations_of_ends > 0
        prev_ent = calculate_entropy_of_ends()
        for i in 1:number_of_iterations_of_ends
            println(i)
            Threads.@threads for var_name in variables_at_ends
                variable_to_factor_messages(variables[var_name])
            end
            Threads.@threads for fact_name in factors_at_ends
                factor_to_variable_messages(factors[fact_name])
            end
            update_all_entropies(variables, variables_at_ends)
            push!(tot_entropy_over_time, total_entropy_of_graph(variables))
            println(tot_entropy_over_time[end])
            cur_ent = calculate_entropy_of_ends()
            if tot_entropy_over_time[end] < 1e-6 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-6 || abs(cur_ent - prev_ent) <= 1e-7
                break
            end
            prev_ent = cur_ent
        end
    end
    read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits, 1) for i in 1:8]

    write(fileFileIOStream, read_off_key == key ? "1\n" : "0\n")
end
close(fileFileIOStream)