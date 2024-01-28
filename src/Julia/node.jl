struct LabelledArray
    array
    axes_labels
end

mutable struct Variable
    name::String
    neighbours::Vector{Factor}
    incoming_messages::Matrix{Float64}
    index_in_neighbours_neighbour::Vector{Int64}
    previous_entropy::Float64
    variables_two_away::Set{String}
    neighbour_index_to_avoid::Int64
    number_of_bits::Int64

    function Variable(name::String, number_of_bits::Int64)
        return new(name, Vector{Factor}[], ones(0, 1 << number_of_bits), Vector{Int64}[], 100_000., Set{String}(), -1, number_of_bits)
    end
end

mutable struct Factor
    name::String
    neighbours::Vector{Variable}
    data::LabelledArray
    incoming_messages::Vector{Vector{Float64}}
    index_in_neighbours_neighbour::Vector{Int64}

    function Factor(name::String, data::LabelledArray)
        return new(name, Vector{Variable}[], data, Vector{Vector{Float64}}[], Vector{Int64}[])
    end
end
