include("node.jl")
include("messages.jl")

function make_xor_prob_table(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 1 << num_of_bits, 1 << num_of_bits)
    for input_a in 1:(1<<num_of_bits)
        for input_b in 1:(1<<num_of_bits)
            output[input_a, input_b, ((input_a-1)⊻(input_b-1))+1] = 1.0
        end
    end
    return output
end

number_of_bits = 16

tab_a = zeros(1 << number_of_bits)
tab_a[1] = 0.8
tab_a[2] = 0.2

p_a = LabelledArray(
    tab_a, ["a"]
)

tab_b = zeros(1 << number_of_bits)
tab_b[1] = 0.7
tab_b[2] = 0.3
p_b = LabelledArray(
    tab_b, ["b"]
)

tab_c = zeros(1 << number_of_bits)
tab_c[1] = 0.9
tab_c[2] = 0.1
p_c = LabelledArray(
    tab_c, ["c"]
)

# xor_table = make_xor_prob_table(number_of_bits)

p_xor = LabelledArray(xor_table, ["a", "b", "c"])

variables = Dict{String,AbsVariable}(
    "a" => Variable{AbsFactor}("a", number_of_bits),
    "b" => Variable{AbsFactor}("b", number_of_bits),
    "c" => Variable{AbsFactor}("c", number_of_bits)
)

factors = Dict{String,AbsFactor}(
    "f1" => Factor{AbsVariable}("f1", p_a),
    "f2" => Factor{AbsVariable}("f2", p_b),
    "f3" => Factor{AbsVariable}("f3", p_c),
    "f_xor" => XorFactor{AbsVariable}("f_xor")
    # "f_xor" => Factor{AbsVariable}("f_xor", p_xor)
)

add_edge_between(variables["a"], factors["f1"])
add_edge_between(variables["a"], factors["f_xor"])

add_edge_between(variables["b"], factors["f2"])
add_edge_between(variables["b"], factors["f_xor"])
add_edge_between(variables["c"], factors["f3"])
add_edge_between(variables["c"], factors["f_xor"])

# The basic idea is to repeated go through and run factor -> variable and variable -> factor messages
for i in 1:4
    for (i, j) in variables
        variable_to_factor_messages(j)
    end

    for (i, j) in factors
        factor_to_variable_messages(j)
    end
end

a_dist = marginal(variables["a"])