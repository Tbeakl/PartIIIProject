struct Variable
    name::String
    neighbours::Array{Factor}
    incoming_messages
end

struct Factor
    name::String
    neighbours::Array{Variable}
    data
    incoming_messages
end

struct LabelledArray
    array
    axes_labels
end