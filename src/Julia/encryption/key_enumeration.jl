using DataStructures
include("rank_estimation.jl")

function make_matrices_of_values(likelihood_tables)
    sorting_permutation = sortperm.(likelihood_tables, rev=true)
    values = [collect(0:(length(likelihood_tables[1]) - 1))[sorting_permutation[i]] for i in eachindex(likelihood_tables)]
    sorted_likelihood_tables = [likelihood_tables[i][sorting_permutation[i]] for i in eachindex(likelihood_tables)]

    return (reduce(hcat, sorted_likelihood_tables), reduce(hcat, values))
end

function is_potential_place_valid(all_explored_areas::Set, potential_location::Vector{Int64}, max_value::Int64)
    for i in eachindex(potential_location)
        if potential_location[i] > max_value
            return false
        elseif potential_location[i] - 1 > 0
            potential_location[i] -= 1
            if !(potential_location in all_explored_areas)
                potential_location[i] += 1
                return false
            end
            potential_location[i] += 1
        end
    end
    return true
end

function get_total_likelihood_of_position(position::Vector{Int64}, sorted_likelihood_matrix::Matrix{Float64})
    likelihood = 0
    @inbounds for i in eachindex(position)
        likelihood += sorted_likelihood_matrix[position[i], i]
    end
    return likelihood
end

function get_position_key_correct(position::Vector{Int64}, sorted_values::Matrix{Int64}, key::Vector{Int64})
    @inbounds for i in eachindex(position)
        if sorted_values[position[i], i] != key[i]
            return false
        end
    end
    return true
end

function key_enumerate(sorted_likelihood_matrix::Matrix{Float64}, sorted_values::Matrix{Int64}, key::Vector{Int64}, iteration_limit::Int64)
    # Need to initialise the priority queue with the initial location and its score
    position = ones(Int64, size(sorted_likelihood_matrix)[2])
    priority_queue = PriorityQueue{Vector, Float64}(Base.Order.Reverse)
    priority_queue[position] = sum(sorted_likelihood_matrix[position, :][1,:])
    all_explored_positions = Set()
    max_value = size(sorted_likelihood_matrix)[1]
    for i in 1:iteration_limit
        if i % 10_000 == 0
            println(i)
            println(peek(priority_queue)[2])
        end
        # println(peek(priority_queue)[2])
        # println("##########")
        # println(all_explored_positions)
        # println("##########")
        position = dequeue!(priority_queue)
        # println(position)
        if get_position_key_correct(position, sorted_values, key)
            return i
        end
        push!(all_explored_positions, position)
        # Need to try increasing the position by 1 in all directions and check if they are valid and if so
        # then add them into the priority queue for further checking
        for i in eachindex(position)
            position[i] += 1
            if is_potential_place_valid(all_explored_positions, position, max_value) && !(position in all_explored_positions)
                # println(position)
                # println(test)
                priority_queue[copy(position)] = get_total_likelihood_of_position(position, sorted_likelihood_matrix)
            end
            position[i] -= 1
        end
    end
    return iteration_limit + 1
end

function turn_key_into_cluster_values(key::Vector{UInt32}, number_of_bits_per_cluster::Int64)
    number_of_clusters = Int64(ceil(32 / number_of_bits_per_cluster))
    values = zeros(Int64, 8 * number_of_clusters)
    for word_number in eachindex(key)
        for i in 1:number_of_clusters
            values[(word_number - 1) * number_of_clusters + i] = (key[word_number] & (((1 << number_of_bits_per_cluster) - 1)) << (number_of_bits_per_cluster * (i - 1))) >> (number_of_bits_per_cluster * (i - 1))
        end
    end
    return values
end

# likelihood_tables = make_log_likelihood_tables_for_key(variables, number_of_bits)
# sorted_likelihood_matrix, sorted_values = make_matrices_of_values(likelihood_tables)
# key_by_cluster = turn_key_into_cluster_values(key, 8)
# calculate_log_likelihood_of_key(key, likelihood_tables, 8)
# actual_rank = key_enumerate(sorted_likelihood_matrix, sorted_values, key_by_cluster, 1 << 24)