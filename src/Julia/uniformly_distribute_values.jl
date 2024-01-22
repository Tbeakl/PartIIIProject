include("messages.jl")

function add_uniform_dist_of_vars(variables::Dict{String, Variable},
    factors::Dict{String, Factor},
    bits_per_cluster::Int64
    )
    for (i,j) in variables
        if occursin("carry", i)
            factor_name = string("f_", i, "_dist_uniform")
            factors[factor_name] = Factor(factor_name, LabelledArray([
                .5
                .5
            ], [i]))
            add_edge_between(j, factors[factor_name])
        elseif occursin("add", i)
            factor_name = string("f_", i, "_dist_uniform")
            factors[factor_name] = Factor(factor_name, LabelledArray(ones(1 << (bits_per_cluster + 1)) ./ (1 << (bits_per_cluster + 1)), [i]))
            add_edge_between(j, factors[factor_name])
        else
            factor_name = string("f_", i, "_dist_uniform")
            factors[factor_name] = Factor(factor_name, LabelledArray(ones(1 << bits_per_cluster) ./ (1 << bits_per_cluster), [i]))
            add_edge_between(j, factors[factor_name])
        end
    end
end

bits_per_cluster = 2
variables = Dict{String, Variable}()
factors = Dict{String, Factor}()
chacha_factor_graph!(variables, factors, bits_per_cluster)
add_uniform_dist_of_vars(variables, factors, bits_per_cluster)

for (i,j) in factors
    factor_to_variable_messages(j)
end

for exec_count in 1:32
    println(exec_count)
    for (i,j) in variables
        variable_to_factor_messages(j)
    end
    for (i,j) in factors
        factor_to_variable_messages(j)
    end
end
