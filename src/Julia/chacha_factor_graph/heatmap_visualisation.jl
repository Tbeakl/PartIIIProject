include("../belief_propagation/node.jl")

function add_in_round(a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64})
    location_execution_counts[a] += 1
    location_execution_counts[c] += 1
end

function xor_in_round(a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64})
    location_execution_counts[d] += 1
    location_execution_counts[b] += 1
end

function first_rot_in_round(a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64},
    bits_per_cluster::Int64)
    if 16 % bits_per_cluster != 0
        location_execution_counts[d] += 1
    end
    if 12 % bits_per_cluster != 0
        location_execution_counts[b] += 1
    end
end

function second_rot_in_round(a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    location_execution_counts::Vector{Int64},
    bits_per_cluster::Int64)
    if 8 % bits_per_cluster != 0
        location_execution_counts[d] += 1
    end
    if 7 % bits_per_cluster != 0
        location_execution_counts[b] += 1
    end
end

function fill_column_index_with_current_values(locations_to_variable_names::Matrix{String},
    location_execution_counts::Vector{Int64},
    clusters_per_word::Int64,
    column_index::Int64,
    run_number::Int64)
    for word_num in 1:16
        for cluster_num in 1:clusters_per_word
            locations_to_variable_names[(word_num - 1) * clusters_per_word + cluster_num, column_index] = string(word_num, "_", location_execution_counts[word_num], "_", cluster_num, "_", run_number)
        end
    end
end

function add_odd_round(locations_to_variable_names::Matrix{String},
    location_execution_counts::Vector{Int64},
    bits_per_cluster::Int64,
    clusters_per_word::Int64,
    initial_column_index::Int64,
    run_number::Int64)
        add_in_round(1, 5, 9, 13, location_execution_counts)
        add_in_round(2, 6, 10, 14, location_execution_counts)
        add_in_round(3, 7, 11, 15, location_execution_counts)
        add_in_round(4, 8, 12, 16, location_execution_counts)

        # Now print to the column_index
        fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index, run_number)

        xor_in_round(1, 5, 9, 13, location_execution_counts)
        xor_in_round(2, 6, 10, 14, location_execution_counts)
        xor_in_round(3, 7, 11, 15, location_execution_counts)
        xor_in_round(4, 8, 12, 16, location_execution_counts)

        # Now print to the column_index
        fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 1, run_number)

        first_rot_in_round(1, 5, 9, 13, location_execution_counts, bits_per_cluster)
        first_rot_in_round(2, 6, 10, 14, location_execution_counts, bits_per_cluster)
        first_rot_in_round(3, 7, 11, 15, location_execution_counts, bits_per_cluster)
        first_rot_in_round(4, 8, 12, 16, location_execution_counts, bits_per_cluster)

        # Now print to the column_index
        fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 2, run_number)

        add_in_round(1, 5, 9, 13, location_execution_counts)
        add_in_round(2, 6, 10, 14, location_execution_counts)
        add_in_round(3, 7, 11, 15, location_execution_counts)
        add_in_round(4, 8, 12, 16, location_execution_counts)

        # Now print to the column_index
        fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 3, run_number)

        xor_in_round(1, 5, 9, 13, location_execution_counts)
        xor_in_round(2, 6, 10, 14, location_execution_counts)
        xor_in_round(3, 7, 11, 15, location_execution_counts)
        xor_in_round(4, 8, 12, 16, location_execution_counts)

        # Now print to the column_index
        fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 4, run_number)

        second_rot_in_round(1, 5, 9, 13, location_execution_counts, bits_per_cluster)
        second_rot_in_round(2, 6, 10, 14, location_execution_counts, bits_per_cluster)
        second_rot_in_round(3, 7, 11, 15, location_execution_counts, bits_per_cluster)
        second_rot_in_round(4, 8, 12, 16, location_execution_counts, bits_per_cluster)

        # Now print to the column_index
        fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 5, run_number)
end

function add_even_round(locations_to_variable_names::Matrix{String},
    location_execution_counts::Vector{Int64},
    bits_per_cluster::Int64,
    clusters_per_word::Int64,
    initial_column_index::Int64,
    run_number::Int64)
    
    add_in_round(1, 6, 11, 16, location_execution_counts)
    add_in_round(2, 7, 12, 13, location_execution_counts)
    add_in_round(3, 8, 9, 14, location_execution_counts)
    add_in_round(4, 5, 10, 15, location_execution_counts)

    # Now print to the column_index
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index, run_number)

    xor_in_round(1, 6, 11, 16, location_execution_counts)
    xor_in_round(2, 7, 12, 13, location_execution_counts)
    xor_in_round(3, 8, 9, 14, location_execution_counts)
    xor_in_round(4, 5, 10, 15, location_execution_counts)

    # Now print to the column_index
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 1, run_number)

    first_rot_in_round(1, 6, 11, 16, location_execution_counts, bits_per_cluster)
    first_rot_in_round(2, 7, 12, 13, location_execution_counts, bits_per_cluster)
    first_rot_in_round(3, 8, 9, 14, location_execution_counts, bits_per_cluster)
    first_rot_in_round(4, 5, 10, 15, location_execution_counts, bits_per_cluster)

    # Now print to the column_index
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 2, run_number)

    add_in_round(1, 6, 11, 16, location_execution_counts)
    add_in_round(2, 7, 12, 13, location_execution_counts)
    add_in_round(3, 8, 9, 14, location_execution_counts)
    add_in_round(4, 5, 10, 15, location_execution_counts)

    # Now print to the column_index
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 3, run_number)

    xor_in_round(1, 6, 11, 16, location_execution_counts)
    xor_in_round(2, 7, 12, 13, location_execution_counts)
    xor_in_round(3, 8, 9, 14, location_execution_counts)
    xor_in_round(4, 5, 10, 15, location_execution_counts)

    # Now print to the column_index
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 4, run_number)

    second_rot_in_round(1, 6, 11, 16, location_execution_counts, bits_per_cluster)
    second_rot_in_round(2, 7, 12, 13, location_execution_counts, bits_per_cluster)
    second_rot_in_round(3, 8, 9, 14, location_execution_counts, bits_per_cluster)
    second_rot_in_round(4, 5, 10, 15, location_execution_counts, bits_per_cluster)

    # Now print to the column_index
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, initial_column_index + 5, run_number)
end

function make_positions_to_var_names(number_of_bits_per_cluster::Int64, run_number::Int64)
    clusters_per_word = Int64(ceil(32 / number_of_bits_per_cluster))
    locations_to_variable_names = fill("", clusters_per_word * 16, 122)
    location_execution_counts = zeros(Int64, 16)
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, 1, run_number)
    for i in 1:10
        add_odd_round(locations_to_variable_names, location_execution_counts, number_of_bits_per_cluster, clusters_per_word, 2 + 12 * (i - 1), run_number)
        add_even_round(locations_to_variable_names, location_execution_counts, number_of_bits_per_cluster, clusters_per_word, 8 + 12 * (i - 1), run_number)
    end
    location_execution_counts[:] .+= 1
    # This is for after the final adds
    fill_column_index_with_current_values(locations_to_variable_names, location_execution_counts, clusters_per_word, 122, run_number)

    ys::Vector{String} = []
    for i in 1:16
        for j in 1:clusters_per_word
            push!(ys, string(i, "_#_", j))
        end
    end

    xs::Vector{String} = ["Initial"]
    for i in 1:20
        push!(xs, string("First Add Round ", i))
        push!(xs, string("First XOR Round ", i))
        push!(xs, string("First Rot Round ", i))
        push!(xs, string("Second Add Round ", i))
        push!(xs, string("Second XOR Round ", i))
        push!(xs, string("Second Rot Round ", i))
    end
    push!(xs, "Final adds")

    return locations_to_variable_names, xs, ys
end

function plot_current_entropy(variables::Dict{String, AbsVariable})
    return function plotting_func(var_name::String)
        return variables[var_name].current_entropy
    end
end

function plot_change_in_entropy(variables::Dict{String, AbsVariable})
    return function plotting_func(var_name::String)
        return abs(variables[var_name].current_entropy - variables[var_name].previous_entropy)
    end
end

function variables_to_heatmap_matrix(variable_names_for_positions::Matrix{String},
    variable_name_to_value_plotting)
    return variable_name_to_value_plotting.(variable_names_for_positions)
end
