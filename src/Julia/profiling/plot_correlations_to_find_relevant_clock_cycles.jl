using HDF5, Plots, Base.Threads, StatsBase, Statistics
plotly()
clock_cycle_sample_number = 405
number_of_intermediate_values = 700
number_of_clock_cycles = (749500 ÷ 500)

all_intermediate_values = zeros(UInt32, length(trace_range_per_file) * length(file_range), number_of_intermediate_values)
mean_value_per_cycle = zeros(Float64, length(trace_range_per_file) * length(file_range), number_of_clock_cycles)

number_of_bits = 8
number_of_clusters = 32 ÷ number_of_bits

data_fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/detection_matrix.hdf5", "r")
all_intermediate_values = read(data_fid["intermediate_values"])
mean_value_per_cycle = read(data_fid["downsampled_matrix"])
close(data_fid)


function calculate_NICV(traces, variances, all_intermediate_values, intermediate_index, number_of_bits, cluster_number)
    intermediate_value_vector = (all_intermediate_values[:, intermediate_index] .>> (number_of_bits * (cluster_number - 1))) .& ((1 << number_of_bits) - 1)
    all_mean_values = zeros(size(traces)[2], 1 << number_of_bits)
    for j in 0:((1 << number_of_bits) - 1)
        all_mean_values[:, j + 1] = mean(traces[intermediate_value_vector .== j, :], dims=1)[1, :]
    end
    return var(all_mean_values, dims=2) ./ variances
end

mean_variances = var(mean_value_per_cycle, dims=1)[1, :]

all_NICV_mean_value = zeros(number_of_intermediate_values * number_of_clusters, size(mean_value_per_cycle)[2])

for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:number_of_clusters
        println(intermediate_value_index, " ", cluster_num)
        intermediate_value_vector = (all_intermediate_values[:, intermediate_value_index] .>> (number_of_bits * (cluster_num - 1))) .& ((1 << number_of_bits) - 1)
        base_prediction_matrix = permutedims(hcat(digits.(intermediate_value_vector, base=2, pad=number_of_bits + 1)...))
        base_prediction_matrix[:, end] .= 1
        for i in 1:(size(mean_value_per_cycle)[2])
            # Need to make the correct weights matrix for the prediction
            β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * mean_value_per_cycle[:, i]
            predicted_values = base_prediction_matrix * β
            all_NICV_mean_value[(number_of_clusters * (intermediate_value_index - 1)) + cluster_num, i] = cor(predicted_values, mean_value_per_cycle[:, i])^2
        end
        # all_NICV_mean_value[(number_of_clusters * (intermediate_value_index - 1)) + cluster_num, :] = calculate_NICV(mean_value_per_cycle, mean_variances, all_intermediate_values, intermediate_value_index, number_of_bits, cluster_num)
    end
end

base_inter_value = 484
p = plot(all_NICV_mean_value[4 * base_inter_value - 3, :])
plot!(p, all_NICV_mean_value[4 * base_inter_value - 2, :])
plot!(p, all_NICV_mean_value[4 * base_inter_value - 1, :])
plot!(p, all_NICV_mean_value[4 * base_inter_value, :])
plot!(p, sum(all_NICV_mean_value[(4 * base_inter_value - 3):(4 * base_inter_value), :], dims=1)[1,:])
hline!(p, [0.005])

# p = plot(mean(Array(all_NICV_mean_value)[8033:8064, :], dims=1)[1, :], legend=false)
# p = plot(Array(all_NICV_mean_value)[8033:8064, :]', legend=false)
# fid = h5open("D:\\ChaChaData\\attack_profiling\\clock_cycles_bitmasks.hdf5", "r")
# cycle_bitmasks = zeros(4, size(max_value_per_cycle)[2])
# for i in 1:4
#     cycle_bitmasks[i, :] = read(fid[string("bitmask_", i + 1004)])
# end
# close(fid)
# plot!(p, cycle_bitmasks' ./ 10, color=:red)

fid = h5open("D:/Year_4_Part_3/Dissertation/SourceCode/PartIIIProject/data/attack_profiling/COR_aligned_traces_2.hdf5", "w")
for intermediate_value_index in 1:number_of_intermediate_values
    for cluster_num in 1:number_of_clusters
        println(intermediate_value_index, " ", cluster_num)
        fid[string("mean_", intermediate_value_index, "_", cluster_num)] = all_NICV_mean_value[(number_of_clusters * (intermediate_value_index - 1)) + cluster_num, :]
    end
end
close(fid)

# fid = h5open("D:\\ChaChaData\\attack_profiling\\downsampled_traces_for_clock_cycle.hdf5", "w")
# fid["mean_values"] = Array(mean_value_per_cycle)
# fid["min_values"] = Array(min_value_per_cycle)
# fid["max_values"] = Array(max_value_per_cycle)
# fid["intermediate_values"] = Array(all_intermediate_values)
# close(fid)