using CryptoSideChannel
using Distributions
include("chacha.jl")

σ = 0.1

function mean_vector_for_value(value)
    output = zeros(Float64, 4)
    output[1] = Base.count_ones(value & 0xFF)
    output[2] = Base.count_ones((value >> 8) & 0xFF)
    output[3] = Base.count_ones((value >> 16) & 0xFF)
    output[4] = Base.count_ones((value >> 24) & 0xFF)
    return output
end


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


# How the full trace (1616) is broken down

# The first 1600 are from the 20 rounds, so that is 80 per round and 20 per quarter round

# Inside a quarter round is broken into 4 parts
# The add has a single part in the trace
# The exclusive-or also has a single part in the trace
# The rotation has three in the trace but I think that only the final should be used
# because that is more fundamental to the actual algorithm

# Last 16 elements from the addition afterwards