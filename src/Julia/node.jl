struct LabelledArray
    array
    axes_labels
end

mutable struct Variable
    name::String
    neighbours::Array{Factor}
    incoming_messages::Dict{String, Union{Float64, Vector{Float64}}}
    previous_entropy::Float64
    variables_two_away::Set{String}
    
    function Variable(name::String)
        return new(name, Array{Factor}[], Dict(), 100_000., Set{String}())
    end
end

mutable struct Factor
    name::String
    neighbours::Array{Variable}
    data::LabelledArray
    incoming_messages::Dict{String, Union{Float64, Vector{Float64}}}
    function Factor(name::String, data::LabelledArray)
        return new(name, Array{Variable}[], data, Dict())
    end
end
