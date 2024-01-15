mutable struct Variable
    name::String
    neighbours::Array{Factor}
    incoming_messages
    
    function Variable(name::String)
        return new(name, Array{Factor}[], Dict())
    end
end

mutable struct Factor
    name::String
    neighbours::Array{Variable}
    data::LabelledArray
    incoming_messages
    function Factor(name::String, data::LabelledArray)
        return new(name, Array{Variable}[], data, Dict())
    end
end

struct LabelledArray
    array
    axes_labels
end