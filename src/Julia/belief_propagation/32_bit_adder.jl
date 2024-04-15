include("node.jl")
include("messages.jl")

function make_add_including_carry_prob_array(num_of_bits::Int64)
    output = zeros(2, 1 << num_of_bits, 1 << num_of_bits, 1 << (num_of_bits + 1))
    for carry_in in 1:2
        for input_a in 1:(1 << num_of_bits)
            for input_b in 1:(1 << num_of_bits)
                output[carry_in, input_a, input_b, input_a + input_b + carry_in - 2] = 1.
            end
        end
    end
    return output
end

function take_top_bit_prob_array(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 2)
    for val_in in 1:(1 << num_of_bits)
        output[val_in, (val_in > 1 << (num_of_bits - 1)) + 1] = 1.
    end
    return output
end

function take_bottom_bits_prob_array(num_of_bits::Int64)
    output = zeros(1 << num_of_bits, 1 << (num_of_bits - 1))
    for val_in in 1:1 << num_of_bits
        output[val_in, ((val_in - 1) % ((1 << (num_of_bits - 1)))) + 1] = 1.
    end
    return output
end

function make_32_bit_adder(number_of_bits_per_cluster::Int64, variables, factors)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    full_add_dist = make_add_including_carry_prob_array(number_of_bits_per_cluster)
    full_add_to_output_dist = take_bottom_bits_prob_array(number_of_bits_per_cluster + 1)
    full_add_carry_dist = take_top_bit_prob_array(number_of_bits_per_cluster + 1)

    variables["carry_0"] = Variable{AbsFactor}("carry_0", 1)
    for i in 1:number_of_clusters
        variables[string("carry_", i)] = Variable{AbsFactor}(string("carry_", i), 1)
        variables[string("input_a_", i)] = Variable{AbsFactor}(string("input_a_", i), number_of_bits_per_cluster)
        variables[string("input_b_", i)] = Variable{AbsFactor}(string("input_b_", i), number_of_bits_per_cluster)
        variables[string("output_temp_", i)] = Variable{AbsFactor}(string("output_temp_", i), number_of_bits_per_cluster + 1)
        variables[string("output_", i)] = Variable{AbsFactor}(string("output_", i), number_of_bits_per_cluster)

        factors[string("f_add_", i)] = Factor{AbsVariable}(string("f_add_", i), LabelledArray(full_add_dist, [string("carry_", i-1), string("input_a_", i), string("input_b_", i), string("output_temp_", i)]))
        factors[string("f_add_output_", i)] = Factor{AbsVariable}(string("f_add_output_", i), LabelledArray(full_add_to_output_dist, [string("output_temp_", i), string("output_", i)]))
        factors[string("f_add_carry_", i)] = Factor{AbsVariable}(string("f_add_carry_", i), LabelledArray(full_add_carry_dist, [string("output_temp_", i), string("carry_", i)]))
    end

    for i in 1:number_of_clusters
        add_edge_between(variables[string("carry_", i - 1)], factors[string("f_add_", i)])
        add_edge_between(variables[string("input_a_", i)], factors[string("f_add_", i)])
        add_edge_between(variables[string("input_b_", i)], factors[string("f_add_", i)])
        add_edge_between(variables[string("output_temp_", i)], factors[string("f_add_", i)])

        add_edge_between(variables[string("output_temp_", i)], factors[string("f_add_output_", i)])
        add_edge_between(variables[string("output_", i)], factors[string("f_add_output_", i)])

        add_edge_between(variables[string("output_temp_", i)], factors[string("f_add_carry_", i)])
        add_edge_between(variables[string("carry_", i)], factors[string("f_add_carry_", i)])
    end
end

function add_base_probabilities(number_of_bits_per_cluster, variables, factors)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    factors["f_carry_0"] = Factor{AbsVariable}("f_carry_0", LabelledArray(
        [
            1.
            0.
        ], ["carry_0"]
    ))
    add_edge_between(variables["carry_0"], factors["f_carry_0"])

    for i in 1:number_of_clusters
        p_input_array = zeros(1 << number_of_bits_per_cluster)
        if i == 1
            p_input_array[2] = .9
            p_input_array[3] = .1
        else
            p_input_array[1] = 1.
        end
        factors[string("f_input_a_", i)] = Factor{AbsVariable}(string("f_input_a_", i), LabelledArray(p_input_array, [string("input_a_", i)]))
        add_edge_between(variables[string("input_a_", i)], factors[string("f_input_a_", i)])

        p_input_array = zeros(1 << number_of_bits_per_cluster)
        if i == 1
            p_input_array[2] = .9
            p_input_array[1] = .1
        else
            p_input_array[1] = 1.
        end
        factors[string("f_input_b_", i)] = Factor{AbsVariable}(string("f_input_b_", i), LabelledArray(p_input_array, [string("input_b_", i)]))
        add_edge_between(variables[string("input_b_", i)], factors[string("f_input_b_", i)])

        p_input_array = zeros(1 << number_of_bits_per_cluster)
        if i == 1
            p_input_array[3] = .9
            p_input_array[4] = .1
        else
            p_input_array[1] = 1.
        end
        factors[string("f_output_", i)] = Factor{AbsVariable}(string("f_output_", i), LabelledArray(p_input_array, [string("output_", i)]))
        add_edge_between(variables[string("output_", i)], factors[string("f_output_", i)])
    end
end

variables = Dict{String, AbsVariable}()
factors = Dict{String, AbsFactor}()

number_of_bits = 2

make_32_bit_adder(number_of_bits, variables, factors)
add_base_probabilities(number_of_bits, variables, factors)

for (i,j) in factors
    factor_to_variable_messages(j)
end

for (i,j) in variables
    variable_to_factor_messages(j)
end

for i in 1:Int64(ceil(32 / number_of_bits))
    println("Output ",i, ": ", marginal(variables[string("output_", i)]))
end