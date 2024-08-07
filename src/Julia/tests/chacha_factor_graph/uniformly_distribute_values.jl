include("../../belief_propagation/node.jl")
include("../../belief_propagation/messages.jl")
include("../../chacha_factor_graph/chacha_factor_graph.jl")

function add_uniform_dist_of_vars(variables::Dict{String, AbsVariable},
    factors::Dict{String, AbsFactor},
    bits_per_cluster::Int64
    )
    for (i,j) in variables
        if occursin("carry", i)
            factor_name = string("f_", i, "_dist_uniform")
            factors[factor_name] = Factor{AbsVariable}(factor_name, LabelledArray([
                .5
                .5
            ], [i]))
            add_edge_between(j, factors[factor_name])
        elseif occursin("add", i)
            factor_name = string("f_", i, "_dist_uniform")
            factors[factor_name] = Factor{AbsVariable}(factor_name, LabelledArray(ones(1 << (bits_per_cluster + 1)) ./ (1 << (bits_per_cluster + 1)), [i]))
            add_edge_between(j, factors[factor_name])
        else
            factor_name = string("f_", i, "_dist_uniform")
            factors[factor_name] = Factor{AbsVariable}(factor_name, LabelledArray(ones(1 << bits_per_cluster) ./ (1 << bits_per_cluster), [i]))
            add_edge_between(j, factors[factor_name])
        end
    end
end

bits_per_cluster = 2
variables = Dict{String, AbsVariable}()
factors = Dict{String, AbsFactor}()
variables_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
factors_by_round::Vector{Set{String}} = [Set{String}() for _ in 1:21]
adds_by_round::Vector{Set{Int64}} = [Set{Int64}() for _ in 1:21]
location_execution_counts = zeros(Int64, 16)
chacha_factor_graph!(variables, factors, bits_per_cluster, variables_by_round, factors_by_round, adds_by_round, location_execution_counts, 1)
add_uniform_dist_of_vars(variables, factors, bits_per_cluster)

all_variables = [keys(variables)...]

for (i,j) in factors
    factor_to_variable_messages(j)
end
update_all_entropies(variables, all_variables)
println("Initial entropy ", total_entropy_of_graph(variables))
for exec_count in 1:30
    for (i,j) in variables
        variable_to_factor_messages(j)
    end
    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    update_all_entropies(variables, all_variables)
    println("Entropy after ", exec_count, " passes ", total_entropy_of_graph(variables))
end
