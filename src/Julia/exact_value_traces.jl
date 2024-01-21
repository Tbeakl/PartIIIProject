using CryptoSideChannel
include("chacha.jl")
include("messages.jl")

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

function add_distribution_of_initial_values(variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    number_of_bits_per_cluster::Int64,
    key::Vector{UInt32}, nonce::Vector{UInt32}, counter::UInt32)
    for i in 1:8
        set_variable_to_value(variables, factors, string(i + 4, "_0"), key[i], number_of_bits_per_cluster)
    end
    set_variable_to_value(variables, factors, "13_0", counter, number_of_bits_per_cluster)
    for i in 1:3
        set_variable_to_value(variables, factors, string(i + 13, "_0"), nonce[i], number_of_bits_per_cluster)
    end
end

function add_value_from_position_in_trace(trace::Vector{Any},
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    variable::Int64,
    amount_to_shift_value_by::Int64
    )
    global position_in_trace
    set_variable_to_value(variables, factors, string(variable, "_", location_execution_counts[variable]), ROTL(trace[position_in_trace], amount_to_shift_value_by), bits_per_cluster)
    position_in_trace += 1
    location_execution_counts[variable] += 1
end

function add_qr_trace(trace::Vector{Any},
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64,
    a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64}
    )
    global position_in_trace

    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, a, 0)
    # Need to do the part with shifting 
    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, bits_per_cluster * (16 ÷ bits_per_cluster))
    position_in_trace += 2
    
    if 16 % bits_per_cluster != 0
        add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, 0)
    else 
        position_in_trace += 1
    end

    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, c, 0)
    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, bits_per_cluster * (12 ÷ bits_per_cluster))
    position_in_trace += 2
    if 12 % bits_per_cluster != 0
        add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, 0)
    else 
        position_in_trace += 1
    end
    

    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, a, 0)
    # Need to do the part with shifting 
    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, bits_per_cluster * (8 ÷ bits_per_cluster))
    position_in_trace += 2
    
    if 8 % bits_per_cluster != 0
        add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, 0)
    else 
        position_in_trace += 1
    end

    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, c, 0)
    add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, bits_per_cluster * (7 ÷ bits_per_cluster))
    position_in_trace += 2
    if 7 % bits_per_cluster != 0
        add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, 0)
    else 
        position_in_trace += 1
    end
end

function add_trace_to_factor_graph(trace::Vector{Any},
    variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64
    )
    global position_in_trace
    position_in_trace = 1
    location_execution_counts = ones(Int64, 16)
    for i in 1:10
        add_qr_trace(trace, variables, factors, bits_per_cluster, 1, 5, 9, 13, location_execution_counts)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 2, 6, 10, 14, location_execution_counts)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 3, 7, 11, 15, location_execution_counts)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 4, 8, 12, 16, location_execution_counts)

        add_qr_trace(trace, variables, factors, bits_per_cluster, 1, 6, 11, 16, location_execution_counts)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 2, 7, 12, 13, location_execution_counts)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 3, 8, 9, 14, location_execution_counts)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 4, 5, 10, 15, location_execution_counts)
    end
    for i in 1:16
        add_value_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, i, 0)
    end
end

# How the full trace (1616) is broken down

# The first 1600 are from the 20 rounds, so that is 80 per round and 20 per quarter round

# Inside a quarter round is broken into 4 parts
# The add has a single part in the trace
# The exclusive-or also has a single part in the trace
# The rotation has three in the trace but I think that only the final should be used
# because that is more fundamental to the actual algorithm

# Last 16 elements from the addition afterwards