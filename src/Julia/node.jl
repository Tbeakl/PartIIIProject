struct Variable
    name::String
    neighbours::Array{Factor}
end

struct Factor
    name::String
    neighbours::Array{Factor}
    data
end

struct LabelledArray
    array
    axes_names
end