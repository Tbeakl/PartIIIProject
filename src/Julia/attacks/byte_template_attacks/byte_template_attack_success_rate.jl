using Plots, Base.Threads, Random, NPZ, FileIO, Dates
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")
include("template_attack_traces.jl")

number_of_bits = 8
dimensions = 8
initial_number_of_iterations = 20
number_of_iterations_of_ends = 400
rounds_for_ends = 5

file_to_write_to_name = "key_templates_nibbles.csv"
fileFileIOStream = open(file_to_write_to_name, "a")

variables = Dict{String,Variable{Factor}}()
factors = Dict{String,Factor{Variable}}()
variables_by_round::Vector{Set{String}} = []
factors_by_round::Vector{Set{String}} = []
adds_by_round::Vector{Vector{Int64}} = []
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
add_starting_constant_values(variables, factors, number_of_bits)

internal_factors = [union(factors_by_round[:]...)...]
internal_variables = [union(variables_by_round[:]...)...]

variables_at_ends = [union(variables_by_round[begin:rounds_for_ends]..., variables_by_round[end- rounds_for_ends - 1:end]...)...]
factors_at_ends = [union(factors_by_round[begin:rounds_for_ends]..., factors_by_round[end- rounds_for_ends - 1:end]...)...]

# Potentially could look at just the very start and end for checking for if we have reached convergence because we 
# don't really need to solve the centre of the graph if the start and end are correct and will not really change any more
number_of_end_rounds_to_check = 2
vars_to_check = [union(variables_by_round[begin:number_of_end_rounds_to_check]..., variables_by_round[end- number_of_end_rounds_to_check - 1:end]...)...]
calculate_entropy_of_ends() = sum([variables[var].current_entropy for var in vars_to_check])

base_path_mean_vectors = "D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\ChaCha_Simulation\\templates_ASCON\\templateLDA_O004\\template_KEY\\template_expect_b"
noise = noise_distribution_fixed_standard_dev(1., dimensions)

for i in 1:1
    variables = Dict{String,Variable{Factor}}()
    factors = Dict{String,Factor{Variable}}()
    variables_by_round::Vector{Set{String}} = []
    factors_by_round::Vector{Set{String}} = []
    adds_by_round::Vector{Vector{Int64}} = []
    location_execution_counts = zeros(Int64, 16)
    chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts)
    add_starting_constant_values(variables, factors, number_of_bits)

    key = generate_random_key()
    nonce = generate_random_nonce()
    counter = generate_random_counter()

    encryption_trace = encrypt_collect_trace(key, nonce, counter, byte_values_for_input)
    encryption_output = encrypt(key, nonce, counter)

    mean_vectors = transpose(npzread(string(base_path_mean_vectors, lpad(string(i % 16), 3, "0"), ".npy")))
    add_byte_template_to_variable = byte_template_value_to_function(mean_vectors, noise)

    # add_values_of_initial_nonce_and_counter(variables, factors, number_of_bits, nonce, counter)
    add_initial_key_dist(variables, factors, number_of_bits, byte_values_for_input.(key), add_byte_template_to_variable)
    add_initial_nonce_and_counter_dist(variables, factors, number_of_bits, byte_values_for_input.(nonce), byte_values_for_input(counter), add_byte_template_to_variable)
    for i in 1:16
        add_byte_template_to_variable(byte_values_for_input(encryption_output[i]), variables, factors, number_of_bits, string(i, "_", location_execution_counts[i]))
        # set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits)
    end

    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits, add_byte_template_to_variable)
    println("Added trace")
    all_variables = [keys(variables)...]
    all_factors = [keys(factors)...]
    tot_entropy_over_time::Vector{Float64} = []
    println(now())
    Threads.@threads for fact_name in all_factors
        factor_to_variable_messages(factors[fact_name])
    end
    println(now())
    for i in 0:adds_by_round[end][end]
        belief_propagate_through_add(variables, factors, number_of_bits, i)
    end

    update_all_entropies(variables, all_variables)
    push!(tot_entropy_over_time, total_entropy_of_graph(variables))

    for i in 1:initial_number_of_iterations
        println(i)
        println(now())
        Threads.@threads for var_name in internal_variables
            variable_to_factor_messages(variables[var_name])
        end
        Threads.@threads for fact_name in internal_factors
            factor_to_variable_messages(factors[fact_name])
        end
        update_all_entropies(variables, all_variables)
        push!(tot_entropy_over_time, total_entropy_of_graph(variables))
        println(tot_entropy_over_time[end])
        if tot_entropy_over_time[end] < 1e-3 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-3
            break
        end
    end

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
        if tot_entropy_over_time[end] < 1e-3 || abs(tot_entropy_over_time[end] - tot_entropy_over_time[end-1]) <= 1e-3 || abs(cur_ent - prev_ent) <= 1e-7
            break
        end
        prev_ent = cur_ent
    end

    read_off_key = [read_most_likely_value_from_variable(variables, string(i + 4, "_0"), number_of_bits) for i in 1:8]
    
    write(fileFileIOStream, read_off_key == key ? "1\n" : "0\n")
end
close(fileFileIOStream)