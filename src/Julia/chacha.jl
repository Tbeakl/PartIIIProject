using CryptoSideChannel
using Distributions

function ROTL(a::T, b) where {T}
    return (a << b) | (a >> (32 - b))
end

function QR(block::Vector{T}, a, b, c, d) where {T}
    block[a] += block[b]
    block[d] ⊻= block[a]
    block[d] = ROTL(block[d], 16)

    block[c] += block[d]
    block[b] ⊻= block[c]
    block[b] = ROTL(block[b], 12)

    block[a] += block[b]
    block[d] ⊻= block[a]
    block[d] = ROTL(block[d], 8)

    block[c] += block[d]
    block[b] ⊻= block[c]
    block[b] = ROTL(block[b], 7)
end

function encrypt(key::Vector{T}, nonce::Vector{T}, counter::T) where {T}
    state = zeros(T, 16)
    state[1] = 0x61707865
    state[2] = 0x3320646e
    state[3] = 0x79622d32
    state[4] = 0x6b206574

    state[5:12] = copy(key)
    state[13] = counter
    state[14:16] = copy(nonce)
    print(state)
    initialState = copy(state)

    for i in 1:10
        QR(state, 1, 5, 9, 13)
        QR(state, 2, 6, 10, 14)
        QR(state, 3, 7, 11, 15)
        QR(state, 4, 8, 12, 16)

        QR(state, 1, 6, 11, 16)
        QR(state, 2, 7, 12, 13)
        QR(state, 3, 8, 9, 14)
        QR(state, 4, 5, 10, 15)
    end

    for i in 1:16
        state[i] += initialState[i]
    end
    return state
end

function mean_vector_for_value(value)
    output = zeros(Float64, 4)
    output[1] = Base.count_ones(value & 0xFF)
    output[2] = Base.count_ones((value >> 8) & 0xFF)
    output[3] = Base.count_ones((value >> 16) & 0xFF)
    output[4] = Base.count_ones((value >> 24) & 0xFF)
    return output
end

σ = 0.1

function noise_distribution(x)
    mean = [0.0, 0.0, 0.0, 0.0]
    C = [σ 0.0 0.0 0.0;
        0.0 σ 0.0 0.0;
        0.0 0.0 σ 0.0;
        0.0 0.0 0.0 σ]

    return MvNormal(mean, C)
end

function encrypt_collect_trace(key::Vector{UInt32}, nonce::Vector{UInt32}, counter::UInt32)
    global trace
    trace = []
    closure = () -> trace

    key_logging = map(x -> Logging.StochasticLog(x, closure, mean_vector_for_value, noise_distribution), key)
    nonce_logging = map(x -> Logging.StochasticLog(x, closure, mean_vector_for_value, noise_distribution), nonce)
    counter_logging = Logging.StochasticLog(counter, closure, mean_vector_for_value, noise_distribution)

    encrypt(key_logging, nonce_logging, counter_logging)

    return copy(trace)
end


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

# Need to work out how the full trace (1616) is broken down

# The first 1600 are from the 20 rounds, so that is 80 per round and 20 per quarter round

# Inside a quarter round is broken into 4 parts
# The add has a single part in the trace
# The exclusive-or also has a single part in the trace
# The rotation has three in the trace but I think that only the final should be used
# because that is more fundamental to the actual algorithm



# Last 16 elements from the addition afterwards