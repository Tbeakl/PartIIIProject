using DataStructures

include("../belief_propagation/node.jl")

function normal_operation(pos_1::Int64,
    location_execution_counts::Vector{Int64})
    location_execution_counts[pos_1] += 1
end

function rotate(pos_1::Int64,
    number_of_bits_to_rotate_by::Int64,
    number_of_bits::Int64,
    location_execution_counts::Vector{Int64})
    if number_of_bits_to_rotate_by % number_of_bits != 0
        location_execution_counts[pos_1] += 1
    end
end

function add_location_to_mapping(location_execution_counts::Vector{Int64},
    number_of_clusters::Int64,
    intermediate_location::Vector{Int64},
    position::Int64,
    mapping::DefaultDict{String,String})
    intermediate_location[1] += 1
    for i in 1:number_of_clusters
        mapping[string(position, "_", location_execution_counts[position], "_", i)] = string(intermediate_location[1], "_", i)
    end
end

function quarter_round(location_execution_counts::Vector{Int64},
    a::Int64,
    b::Int64,
    c::Int64,
    d::Int64,
    number_of_bits::Int64,
    number_of_clusters::Int64,
    intermediate_location::Vector{Int64},
    mapping::DefaultDict{String,String})

    normal_operation(a, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, a, mapping)
    normal_operation(d, location_execution_counts)
    rotate(d, 16, number_of_bits, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, d, mapping)

    normal_operation(c, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, c, mapping)
    normal_operation(b, location_execution_counts)
    rotate(b, 12, number_of_bits, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, b, mapping)

    normal_operation(a, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, a, mapping)
    normal_operation(d, location_execution_counts)
    rotate(d, 8, number_of_bits, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, d, mapping)


    normal_operation(c, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, c, mapping)
    normal_operation(b, location_execution_counts)
    rotate(b, 7, number_of_bits, location_execution_counts)
    add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, b, mapping)
end

function turn_intermediate_name_to_intermediate_index(number_of_bits::Int64)
    mapping = DefaultDict{String, String}("NAN")
    location_execution_counts::Vector{Int64} = zeros(Int64, 16)

    number_of_clusters::Int64 = 32 ÷ number_of_bits

    for intermediate_value in 1:8
        for cluster_num in 1:number_of_clusters
            mapping[string(intermediate_value + 4, "_0_", cluster_num)] = string(intermediate_value, "_", cluster_num)
        end
    end

    for intermediate_value in 9:11
        for cluster_num in 1:number_of_clusters
            mapping[string(intermediate_value + 5, "_0_", cluster_num)] = string(intermediate_value, "_", cluster_num)
        end
    end

    for cluster_num in 1:number_of_clusters
        mapping[string("13_0_", cluster_num)] = string(12, "_", cluster_num)
    end
    intermediate_location::Vector{Int64} = [28]

    for i in 1:10
        quarter_round(location_execution_counts, 1, 5, 9, 13, number_of_bits, number_of_clusters, intermediate_location, mapping)
        quarter_round(location_execution_counts, 2, 6, 10, 14, number_of_bits, number_of_clusters, intermediate_location, mapping)
        quarter_round(location_execution_counts, 3, 7, 11, 15, number_of_bits, number_of_clusters, intermediate_location, mapping)
        quarter_round(location_execution_counts, 4, 8, 12, 16, number_of_bits, number_of_clusters, intermediate_location, mapping)

        quarter_round(location_execution_counts, 1, 6, 11, 16, number_of_bits, number_of_clusters, intermediate_location, mapping)
        quarter_round(location_execution_counts, 2, 7, 12, 13, number_of_bits, number_of_clusters, intermediate_location, mapping)
        quarter_round(location_execution_counts, 3, 8, 9, 14, number_of_bits, number_of_clusters, intermediate_location, mapping)
        quarter_round(location_execution_counts, 4, 5, 10, 15, number_of_bits, number_of_clusters, intermediate_location, mapping)
    end

    for i in 1:16
        normal_operation(i, location_execution_counts)
        add_location_to_mapping(location_execution_counts, number_of_clusters, intermediate_location, i, mapping)
    end
    return mapping
end

function map_to_values(location::String, values::Vector, number_of_templates_per_value::Int64)
    if location == "NAN"
        return Float64(NaN)
    else
        split_location = split(location, "_")
        return Float64(values[number_of_templates_per_value * (parse(Int64, split_location[1]) - 1) + parse(Int64, split_location[2])])
    end
end
