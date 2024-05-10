using HDF5, Plots, Base.Threads, StatsBase, Statistics, CUDA
plotly()
function main()
    clock_cycle_sample_number = 46
    number_of_intermediate_values = 700
    number_of_clock_cycles = (799850 ÷ 50)

    number_of_bits = 8
    number_of_clusters = 32 ÷ number_of_bits

    path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

    data_fid = h5open(string(path_to_data, "attack_profiling/8_on_32/detection.hdf5"), "r")
    all_intermediate_values = read(data_fid["intermediate_values"])
    mean_value_per_cycle = read(data_fid["downsampled_matrix"])
    close(data_fid)

    function calculate_linear_correlation(base_prediction_matrix, values_per_cycle, i)
        β = (base_prediction_matrix' * base_prediction_matrix) \ base_prediction_matrix' * values_per_cycle[:, i]
        predicted_values = base_prediction_matrix * β
        return predicted_values
    end

    normal_memory_mean_val_per_cycle = eachcol(mean_value_per_cycle)
    mean_value_per_cycle = CuArray(mean_value_per_cycle)

    all_NICV_mean_value = zeros(number_of_intermediate_values * number_of_clusters, size(mean_value_per_cycle)[2])


    for intermediate_value_index in 1:number_of_intermediate_values
        for cluster_num in 1:number_of_clusters
            println(intermediate_value_index, " ", cluster_num)
            intermediate_value_vector = (all_intermediate_values[:, intermediate_value_index] .>> (number_of_bits * (cluster_num - 1))) .& ((1 << number_of_bits) - 1)
            base_prediction_matrix = permutedims(hcat(digits.(intermediate_value_vector, base=2, pad=number_of_bits + 1)...))
            base_prediction_matrix[:, end] .= 1
            base_prediction_matrix = CuArray(Int32.(base_prediction_matrix))
            predicted_vectors = calculate_linear_correlation.(Ref(base_prediction_matrix), Ref(mean_value_per_cycle), 1:(size(mean_value_per_cycle)[2]))
            all_NICV_mean_value[(number_of_clusters*(intermediate_value_index-1))+cluster_num, :] = cor.(normal_memory_mean_val_per_cycle, Array.(predicted_vectors)) .^ 2
        end
    end
    # base_inter_value = 251
    # p = plot(all_NICV_mean_value[4 * base_inter_value - 3, :])
    # plot!(p, all_NICV_mean_value[4 * base_inter_value - 2, :])
    # plot!(p, all_NICV_mean_value[4 * base_inter_value - 1, :])
    # plot!(p, all_NICV_mean_value[4 * base_inter_value, :])
    # # plot!(p, mean(all_NICV_mean_value[(4 * base_inter_value - 3):(4 * base_inter_value), :], dims=1)[1,:])
    # hline!(p, [0.004])


    fid = h5open(path_to_data * "attack_profiling/8_on_32/COR_aligned_traces.hdf5", "w")
    for intermediate_value_index in 1:number_of_intermediate_values
        for cluster_num in 1:number_of_clusters
            println(intermediate_value_index, " ", cluster_num)
            fid[string("mean_", intermediate_value_index, "_", cluster_num)] = all_NICV_mean_value[(number_of_clusters*(intermediate_value_index-1))+cluster_num, :]
        end
    end
    close(fid)
end

main()