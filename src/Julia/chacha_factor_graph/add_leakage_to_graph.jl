include("../belief_propagation/node.jl")

# Currently need the clusters to fall exactly along with the leakages
function add_distribution_from_position_in_trace(trace::Vector{Any},
    variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    bits_per_cluster::Int64,
    location_execution_counts::Vector{Int64},
    variable::Int64,
    run_number::Int64,
    trace_value_to_graph_function
    )
    global position_in_trace
    
    trace_value_to_graph_function(trace[position_in_trace], variables, factors, bits_per_cluster, string(variable, "_", location_execution_counts[variable]), run_number)
    # Need to pay attention to how to deal with rotations
    position_in_trace += 1
    location_execution_counts[variable] += 1
end

function add_qr_trace(trace::Vector{Any},
    variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    bits_per_cluster::Int64,
    a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64},
    run_number::Int64,
    trace_value_to_graph_function_add_output,
    trace_value_to_graph_function_rot_output
    )
    global position_in_trace

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, a, run_number, trace_value_to_graph_function_add_output)
    # Need to do the part with shifting
    position_in_trace += 3
    if 16 % bits_per_cluster != 0
        # This means that the result of the rotation is directly going to be the output after the rotation
        location_execution_counts[d] += 1 
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, run_number, trace_value_to_graph_function_rot_output)

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, c, run_number, trace_value_to_graph_function_add_output)
    position_in_trace += 3
    if 12 % bits_per_cluster != 0
        location_execution_counts[b] += 1
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, run_number, trace_value_to_graph_function_rot_output)

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, a, run_number, trace_value_to_graph_function_add_output)
    # Need to do the part with shifting
    position_in_trace += 3
    if 8 % bits_per_cluster != 0
        # This means that the result of the rotation is directly going to be the output after the rotation
        location_execution_counts[d] += 1 
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, d, run_number, trace_value_to_graph_function_rot_output)

    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, c, run_number, trace_value_to_graph_function_add_output)
    position_in_trace += 3
    if 7 % bits_per_cluster != 0
        location_execution_counts[b] += 1
    end
    add_distribution_from_position_in_trace(trace, variables, factors, bits_per_cluster, location_execution_counts, b, run_number, trace_value_to_graph_function_rot_output)
end

function add_trace_to_factor_graph(trace::Vector{Any},
    variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    bits_per_cluster::Int64,
    run_number::Int64,
    trace_value_to_graph_function_add_output,
    trace_value_to_graph_function_rot_output
    )
    global position_in_trace
    position_in_trace = 1
    location_execution_counts = ones(Int64, 16)
    for i in 1:10
        add_qr_trace(trace, variables, factors, bits_per_cluster, 1, 5, 9, 13, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 2, 6, 10, 14, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 3, 7, 11, 15, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 4, 8, 12, 16, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)

        add_qr_trace(trace, variables, factors, bits_per_cluster, 1, 6, 11, 16, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 2, 7, 12, 13, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 3, 8, 9, 14, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
        add_qr_trace(trace, variables, factors, bits_per_cluster, 4, 5, 10, 15, location_execution_counts, run_number, trace_value_to_graph_function_add_output, trace_value_to_graph_function_rot_output)
    end
    # Need to have another part for setting the values to be exactly what they were in the output of the function because
    # we assume that we have access to those actual values compared to just their leakage values
end

function add_initial_key_dist(variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    bits_per_cluster::Int64,
    key_leakage::Vector,
    run_number::Int64,
    add_dist_to_variable)
    for i in 1:8
        add_dist_to_variable(key_leakage[i], variables, factors, bits_per_cluster, string(i + 4, "_", 0), run_number)
    end
end

function add_initial_nonce_and_counter_dist(variables::Dict{String, Variable{Factor}},
    factors::Dict{String, Factor{Variable}},
    number_of_bits_per_cluster::Int64,
    nonce_leakage::Vector, 
    counter_leakage,
    run_number::Int64,
    add_dist_to_variable)
    add_dist_to_variable(counter_leakage, variables, factors, number_of_bits_per_cluster, string("13_0"), run_number)
    for i in 1:3
        add_dist_to_variable(nonce_leakage[i], variables, factors, number_of_bits_per_cluster, string(i + 13, "_0"), run_number)
    end
end