using HDF5, Plots
include("common_functions.jl")
number_of_intermediate_values = 700
number_of_bits = 8
templates_per_intermediate_value = 32 ÷ number_of_bits
number_of_clock_cycles = 1499

all_bitmasks = zeros(UInt8, number_of_clock_cycles, number_of_intermediate_values * templates_per_intermediate_value)

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_no_dilation_cor.hdf5", "r")
for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:templates_per_intermediate_value
        println(intermediate_value_index, " ", cluster_num)
        all_bitmasks[:, (templates_per_intermediate_value*(intermediate_value_index-1))+cluster_num] = read(fid[string("bitmask_", intermediate_value_index, "_", cluster_num)])
    end
end
close(fid)

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/clock_cycles_bitmasks_dilated_by_leakage_event.hdf5", "w")
for intermediate_value in 1:number_of_intermediate_values

    current_bitmap = Int8.(all_bitmasks[:, templates_per_intermediate_value*(intermediate_value-1)+1])
    current_bitmap = dilate_after(current_bitmap, 2)
    current_bitmap = dilate_infront(current_bitmap, 4)
    plot(current_bitmap)

    # Find all the seperate leakage events outputting each one as a seperate bitmask which can then be used
    # to build templates in the future which hopefully when combined give better outcomes than a single template
    # which takes all of the versions into account because hopefully they will have slightly different lots of noise
    # they are subject to, this I am hoping with especially for repeated runs with the same value

    diff_values = diff(current_bitmap)
    starting_indices = findall(x -> x == 1, diff_values) .+ 1
    ending_indices = findall(x -> x == -1, diff_values)
    for i in eachindex(starting_indices)
        leakage_bitmap = zeros(UInt8, number_of_clock_cycles)
        leakage_bitmap[starting_indices[i]:ending_indices[i]] .= 1
        fid[string("bitmask_", intermediate_value, "_", i)] = leakage_bitmap
    end
end
close(fid)