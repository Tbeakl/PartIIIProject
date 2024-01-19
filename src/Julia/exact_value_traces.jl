using CryptoSideChannel
include("chacha.jl")

function encrypt_collect_trace_full_values(key::Vector{UInt32}, nonce::Vector{UInt32}, counter::UInt32)
    global trace
    trace = []
    closure = () -> trace

    key_logging = map(x -> Logging.FullLog(x, closure), key)
    nonce_logging = map(x -> Logging.FullLog(x, closure), nonce)
    counter_logging = Logging.FullLog(counter, closure)

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