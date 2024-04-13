using Plots, Base.Threads, Random, BenchmarkTools
include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")
include("../../chacha_factor_graph/add_leakage_to_graph.jl")
include("../../chacha_factor_graph/heatmap_visualisation.jl")
include("../../encryption/leakage_functions.jl")
include("../../encryption/chacha.jl")

# key = zeros(UInt32, 8)
# nonce = zeros(UInt32, 3)
# counter::UInt32 = 0

key = [0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c]
nonce = [0x09000000, 0x4a000000, 0x00000000]
counter::UInt32 = 1

number_of_bits = 4

# key = generate_random_key()
# nonce = generate_random_nonce()
# counter = generate_random_counter()

encryption_output = encrypt(key, nonce, counter)

variables = Dict{String,Variable{Factor}}()
factors = Dict{String,Factor{Variable}}()
variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, number_of_bits, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)
add_starting_constant_values(variables, factors, number_of_bits, 1)
add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter, 1)

for i in 1:16
    set_variable_to_value(variables, factors, string(i, "_", location_execution_counts[i]), encryption_output[i], number_of_bits, 1)
end

all_variables = [keys(variables)...]
all_factors = [keys(factors)...]

heatmap_plotting_function = plot_current_entropy(variables)
visualisation_of_entropy::Vector{Matrix{Float64}} = []
visualisation_variables, x_labels, y_labels = make_positions_to_var_names(number_of_bits, 1)
tot_entropy_over_time::Vector{Float64} = []

# # Perform the dynamic message passing around the graph to push information through the entire graph
# dynamic_belief_propogate_through_graph(variables, factors, 100_000)

# This seems to heavily converge to bad conditions, may be better to do some more full passes at some point to hopefully get to a correct convergence

# There is some issue with the passing of the probabilities in this because they are all going to zeros
# which is a major problem, not really sure where it is coming from because of the need 

Threads.@threads for fact_name in all_factors
    factor_to_variable_messages(factors[fact_name])
end

# Think it could be helpful to initially push through the values which have actually been set away from uniform
# by the leakage, because there are parts which are left unknown which may not be ideal to push through immediately because
# they are providing no information about what the key actually is - I will currently just do this for adds because
# I think these will be the main source of the issues although there could be ones also for the rotations probably should 
# just add the correct leakages in for before the rotation has taken place as well because that might actually be recoverable

# for add_nums in adds_by_round
#     for add_num in add_nums
#         belief_propagate_through_add(variables, factors, number_of_bits, add_num, 1)
#     end
# end

update_all_entropies(variables, all_variables)
push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))
push!(tot_entropy_over_time, total_entropy_of_graph(variables))
println(tot_entropy_over_time[end])

internal_factors = [union(factors_by_round[:]...)...]
internal_variables = [union(variables_by_round[:]...)...]

initial_number_of_iterations = 250

for i in 1:initial_number_of_iterations
    println(i)
    Threads.@threads for var_name in internal_variables
        variable_to_factor_messages(variables[var_name])
    end
    Threads.@threads for fact_name in internal_factors
        factor_to_variable_messages(factors[fact_name])
    end
    update_all_entropies(variables, all_variables)
    push!(visualisation_of_entropy, variables_to_heatmap_matrix(visualisation_variables, heatmap_plotting_function))
    
    push!(tot_entropy_over_time, total_entropy_of_graph(variables))
    println(tot_entropy_over_time[end])
    if tot_entropy_over_time[end] == 0
        break
    end
end

@benchmark begin
    for var_name in internal_variables
        variable_to_factor_messages(variables[var_name])
    end
    for fact_name in internal_factors
        factor_to_variable_messages(factors[fact_name])
    end
end

anim = @animate for i in eachindex(visualisation_of_entropy)
    heatmap(visualisation_of_entropy[i]; title=string("Round ", i - 1, " entropy of variables"), clim=(0, number_of_bits)) # 
end

# heatmap(visualisation_of_entropy[1]; title=string("Round ", 0, " entropy of variables")) # clim=(0, number_of_bits),
gif(anim, string("test_test_", number_of_bits, ".gif"), fps=10)