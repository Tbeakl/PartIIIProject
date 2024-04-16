struct LabelledArray
    array
    axes_labels
end

abstract type AbsVariable end

abstract type AbsFactor end

mutable struct Variable{T} <: AbsVariable
    name::String
    neighbours::Vector{T}
    incoming_messages::Matrix{Float64}
    index_in_neighbours_neighbour::Vector{Int64}
    current_entropy::Float64
    previous_entropy::Float64
    variables_two_away::Set{String}
    neighbour_index_to_avoid::Int64
    number_of_bits::Int64

    function Variable{T}(name::String, number_of_bits::Int64) where {T}
        return new(name, Vector{T}[], ones(0, 1 << number_of_bits), Vector{Int64}[], number_of_bits, 100_000., Set{String}(), -1, number_of_bits)
    end
end

mutable struct Factor{T} <: AbsFactor
    name::String
    neighbours::Vector{T}
    data::LabelledArray
    incoming_messages::Vector{Vector{Float64}}
    index_in_neighbours_neighbour::Vector{Int64}

    function Factor{T}(name::String, data::LabelledArray) where {T}
        return new(name, Vector{T}[], data, Vector{Vector{Float64}}[], Vector{Int64}[])
    end
end

mutable struct XorFactor{T} <: AbsFactor
    name::String
    neighbours::Vector{T}
    incoming_messages::Vector{Vector{Float64}}
    index_in_neighbours_neighbour::Vector{Int64}

    function XorFactor{T}(name::String) where {T}
        return new(name, Vector{T}[], Vector{Vector{Float64}}[], Vector{Int64}[])
    end
end

mutable struct AddFactor{T} <: AbsFactor
    name::String
    neighbours::Vector{T}
    incoming_messages::Vector{Vector{Float64}}
    index_in_neighbours_neighbour::Vector{Int64}

    function AddFactor{T}(name::String) where {T}
        return new(name, Vector{T}[], Vector{Vector{Float64}}[], Vector{Int64}[])
    end
end