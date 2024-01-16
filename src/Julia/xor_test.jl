include("messages.jl")

p_a = LabelledArray(
    [.8
     .2], ["a"]
)

p_b = LabelledArray(
    [.7
     .3], ["b"]
)

p_c = LabelledArray(
    [
    .9
    .1
    ], ["c"]
)

xor_table = zeros(2,2,2)

xor_table[1, :, :] = [1. 0.
                      0. 1.]

xor_table[2, :, :] = [0. 1.
                      1. 0.]

p_xor = LabelledArray(xor_table, ["a", "b", "c"])

variables = Dict(
    "a" => Variable("a"),
    "b" => Variable("b"),
    "c" => Variable("c")
)

factors = Dict(
    "f1" => Factor("f1", p_a),
    "f2" => Factor("f2", p_b),
    "f3" => Factor("f3", p_c),
    "f_xor" => Factor("f_xor", p_xor)
)

add_edge_between(variables["a"], factors["f1"])
add_edge_between(variables["a"], factors["f_xor"])

add_edge_between(variables["b"], factors["f2"])
add_edge_between(variables["b"], factors["f_xor"])
add_edge_between(variables["c"], factors["f3"])
add_edge_between(variables["c"], factors["f_xor"])

# The basic idea is to repeated go through and run factor -> variable and variable -> factor messages

for (i,j) in factors
    factor_to_variable_messages(j)
end

for (i,j) in variables
    variable_to_factor_messages(j)
end

a_dist = marginal(variables["a"])