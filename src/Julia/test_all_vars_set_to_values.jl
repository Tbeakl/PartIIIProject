using Test
include("chacha_factor_graph.jl")
include("exact_value_traces.jl")

function all_zeros_single_bit_clusters()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 1
    variables = Dict{String, Variable}()
    factors = Dict{String, Factor}()
    chacha_factor_graph!(variables, factors, number_of_bits)
    add_starting_constant_values(variables, factors, number_of_bits)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter)
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits)

    for execution_count in 1:32
        println(execution_count)
        for (i,j) in factors
            factor_to_variable_messages(j)
        end
        for (i,j) in variables
            variable_to_factor_messages(j)
        end
    end
    # Make a copy of all of the distributions to check for a fixed point
    variable_marginal_distributions = Dict{String, Vector{Float64}}()
    for (i,j) in variables
        variable_marginal_distributions[i] = marginal(j)
    end

    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end

    for (i,j) in variables
        @test marginal(j) == variable_marginal_distributions[i]
    end
end

function all_zeros_2_bit_clusters()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 2
    variables = Dict{String, Variable}()
    factors = Dict{String, Factor}()
    chacha_factor_graph!(variables, factors, number_of_bits)
    add_starting_constant_values(variables, factors, number_of_bits)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter)
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits)

    for execution_count in 1:16
        println(execution_count)
        for (i,j) in factors
            factor_to_variable_messages(j)
        end
        for (i,j) in variables
            variable_to_factor_messages(j)
        end
    end
    # Make a copy of all of the distributions to check for a fixed point
    variable_marginal_distributions = Dict{String, Vector{Float64}}()
    for (i,j) in variables
        variable_marginal_distributions[i] = marginal(j)
    end

    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end

    for (i,j) in variables
        @test marginal(j) == variable_marginal_distributions[i]
    end
end

function all_zeros_4_bit_clusters()
    key = zeros(UInt32, 8)
    nonce = zeros(UInt32, 3)
    counter :: UInt32 = 0
    encryption_trace = encrypt_collect_trace_full_values(key, nonce, counter)
    
    number_of_bits = 4
    variables = Dict{String, Variable}()
    factors = Dict{String, Factor}()
    chacha_factor_graph!(variables, factors, number_of_bits)
    add_starting_constant_values(variables, factors, number_of_bits)
    add_distribution_of_initial_values(variables, factors, number_of_bits, key, nonce, counter)
    add_trace_to_factor_graph(encryption_trace, variables, factors, number_of_bits)

    for execution_count in 1:8
        println(execution_count)
        for (i,j) in factors
            factor_to_variable_messages(j)
        end
        for (i,j) in variables
            variable_to_factor_messages(j)
        end
    end
    # Make a copy of all of the distributions to check for a fixed point
    variable_marginal_distributions = Dict{String, Vector{Float64}}()
    for (i,j) in variables
        variable_marginal_distributions[i] = marginal(j)
    end

    for (i,j) in factors
        factor_to_variable_messages(j)
    end
    for (i,j) in variables
        variable_to_factor_messages(j)
    end

    for (i,j) in variables
        @test marginal(j) == variable_marginal_distributions[i]
    end
end

all_zeros_single_bit_clusters()
all_zeros_2_bit_clusters()
all_zeros_4_bit_clusters()