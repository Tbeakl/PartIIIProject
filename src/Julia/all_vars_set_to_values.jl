include("chacha_factor_graph.jl")
include("exact_value_traces.jl")

base_key = zeros(UInt32, 8)
base_nonce = zeros(UInt32, 3)
base_counter::UInt32 = 0

encryption_trace = encrypt_collect_trace_full_values(base_key, base_nonce, base_counter)

number_of_bits = 2
variables = Dict{String, Variable}()
factors = Dict{String, Factor}()
chacha_factor_graph!(variables, factors, number_of_bits)
add_starting_constant_values(variables, factors, number_of_bits)
add_distribution_of_initial_values(variables, factors, number_of_bits, base_key, base_nonce, base_counter)

add_trace_to_factor_graph(trace, variables, factors, number_of_bits)

for execution_number in 1:100
    println(execution_number)
    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end
end


for (i,j) in variables
    variable_to_factor_messages(j)
end
for (i,j) in factors
    factor_to_variable_messages(j)
end

println("1_0: ", read_most_likely_value_from_variable(variables, "1_0", number_of_bits))
println("5_0: ", read_most_likely_value_from_variable(variables, "5_0", number_of_bits))
println("9_0: ", read_most_likely_value_from_variable(variables, "9_0", number_of_bits))
println("13_0: ", read_most_likely_value_from_variable(variables, "13_0", number_of_bits))

println("1_1: ", read_most_likely_value_from_variable(variables, "1_1", number_of_bits))
println("5_1: ", read_most_likely_value_from_variable(variables, "5_1", number_of_bits))
println("9_1: ", read_most_likely_value_from_variable(variables, "9_1", number_of_bits))
println("13_1: ", read_most_likely_value_from_variable(variables, "13_1", number_of_bits))

# println(factors["f_add_1_1"].incoming_messages)
# println(factors["f_xor_1_1"].incoming_messages)
# println(marginal(variables["1_1_1"]))