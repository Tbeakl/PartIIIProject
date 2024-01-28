using DataStructures
using NaNMath

include("node.jl")
include("messages.jl")

calculate_entropy(prob_dist) = -NaNMath.sum(prob_dist .* log2.(prob_dist))

function populate_with_variables_two_away(variables::Dict{String, Variable})
    for (var_name, variable) in variables
        all_vars_two_away = Set{String}()
        for factor in variable.neighbours
            union!(all_vars_two_away, [var.name for var in factor.neighbours])
        end
        variable.variables_two_away = all_vars_two_away
    end
end

function update_variable_priority(variable::Variable,
    priority_queue::PriorityQueue{String, Float64})
    prob_dist = marginal(variable)
    new_entropy = calculate_entropy(prob_dist) 
    if isnan(variable.previous_entropy) && isnan(new_entropy)
        priority_queue[variable.name] = 0.
    elseif isnan(new_entropy)
        priority_queue[variable.name] = Float64(length(prob_dist))
    else
        priority_queue[variable.name] = variable.previous_entropy - new_entropy
    end
end

function variable_to_factor_messages_dynamic_scheduling(variables::Dict{String, Variable},
    priority_queue::PriorityQueue{String, Float64})
    variable_name_to_update = peek(priority_queue)[1]
    variables[variable_name_to_update].previous_entropy = calculate_entropy(marginal(variables[variable_name_to_update]))
    
    variable_to_factor_messages(variables[variable_name_to_update])
    
    factor_to_variable_messages.(variables[variable_name_to_update].neighbours)
    # Update all the priorities to be the current amount of change in entropy since the last update
    for var_name_to_update_change in variables[variable_name_to_update].variables_two_away
        update_variable_priority(variables[var_name_to_update_change], priority_queue)
    end
end

function dynamic_belief_propogate_through_graph(variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    max_number_of_updates::Int64
    )
    populate_with_variables_two_away(variables)
    println("Populated vars two away")
    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    println("Passed initial messages")
    priority_queue = PriorityQueue{String, Float64}(Base.Order.Reverse)
    for (i,j) in variables
        update_variable_priority(j, priority_queue)
    end
    println("Initial entropy in graph ", total_entropy_of_graph(variables))
    println("Filled the priority queue")
    println("Initial first element in priority queue ", peek(priority_queue))
    number_of_updates::Int64 = 0
    while number_of_updates < max_number_of_updates && abs(peek(priority_queue)[2]) > 1e-6
        # println(peek(priority_queue))
        variable_to_factor_messages_dynamic_scheduling(variables, priority_queue)
        # println("Element after update ", peek(priority_queue))
        number_of_updates += 1
        if number_of_updates & 8191 == 0
            println("Entropy remaining in graph after ", number_of_updates, " updates is ", total_entropy_of_graph(variables))
        end
    end
    println("Took ", number_of_updates, " updates")
end
