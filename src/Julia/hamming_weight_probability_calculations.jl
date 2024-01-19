
function probability_of_choosing_without_replacement(number_in_input::Int64, number_set_in_input::Int64, number_chosen::Int64, number_set_in_output::Int64)
    if number_set_in_output > number_set_in_input
        return 0.0
    end
    output_numerator = 1.0
    output_denominator = 1.0
    for i in 1:number_chosen
        if number_set_in_output > 0
            output_numerator *= number_set_in_input
            number_set_in_input -= 1
            number_set_in_output -= 1
        else
            output_numerator *= (number_in_input - number_set_in_input)
        end
        output_denominator *= number_in_input
        number_in_input -= 1
    end
    return output_numerator / output_denominator
end

# Only support the splitting up of the parts when it is split evenly between the different parts
function probability_of_hamming_values(standard_deviation::Float64, number_of_bits_in_input::Int64, number_of_bits_in_output::Int64, leakage_value::Float64)
    number_of_set_bits_in_original = 0:number_of_bits_in_input
    number_of_set_bits_in_output = 0:number_of_bits_in_output

    dist = Normal(leakage_value, standard_deviation)
    likelihoods_of_set_bits = pdf(dist, number_of_set_bits_in_original)

    likelihood_table = zeros(number_of_bits_in_input + 1, number_of_bits_in_output + 1)
    for i in number_of_set_bits_in_output
        likelihood_table[:, i+1] = likelihoods_of_set_bits .* probability_of_choosing_without_replacement.(number_of_bits_in_input, number_of_set_bits_in_original, number_of_bits_in_output, i)
    end
    return sum(likelihood_table, dims=1)
end
