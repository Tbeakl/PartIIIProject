using CryptoSideChannel
using Random
include("chacha.jl")

function byte_hamming_weight_for_value(value)
    output = zeros(Int64, 4)
    output[1] = Base.count_ones(value & 0xFF)
    output[2] = Base.count_ones((value >> 8) & 0xFF)
    output[3] = Base.count_ones((value >> 16) & 0xFF)
    output[4] = Base.count_ones((value >> 24) & 0xFF)
    return output
end

function byte_values_for_input(value)
    output = zeros(Int64, 4)
    output[1] = value & 0xFF
    output[2] = (value >> 8) & 0xFF
    output[3] = (value >> 16) & 0xFF
    output[4] = (value >> 24) & 0xFF
    return output
end

full_value(value) = value

# Think I should probably add on some other factor for determining how the noise should be shaped
# but it might not matter because of the fact that we set the signal to noise ratio when generating the
# mean vectors
function noise_distribution_random_variances(dimensions::Int64)
    rng = MersenneTwister()
    # This is making the noise independent in all directions could experiment
    # with correlating it
    cov = rand(rng, dimensions)
    return MvNormal(cov)
end

function noise_distribution_fixed_standard_dev(standard_deviation)
    variance = standard_deviation * standard_deviation
    covariance_matrix = ones(4) .* variance
    return MvNormal(covariance_matrix)
end

function encrypt_collect_trace(key::Vector{UInt32}, nonce::Vector{UInt32}, counter::UInt32, leakage_function)
    global trace
    trace = []
    closure = () -> trace
    noise_dist_closure(x) = noise_distribution
    key_logging = map(x -> Logging.SingleFunctionLog(x, closure, leakage_function), key)
    nonce_logging = map(x -> Logging.SingleFunctionLog(x, closure, leakage_function), nonce)
    counter_logging = Logging.SingleFunctionLog(counter, closure, leakage_function)

    encrypt(key_logging, nonce_logging, counter_logging)

    return copy(trace)
end
