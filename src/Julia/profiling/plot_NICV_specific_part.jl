using HDF5, MultivariateStats, Plots, StatsBase, Statistics, LinearAlgebra, Images, CUDA
plotly()
fid = h5open("D:\\ChaChaData\\attack_profiling\\original_1002_values.hdf5", "r")
intermediate_value_vector = read(fid["intermediate_values"])
original_matrix_of_current_data = Int16.(fid["downsampled_matrix"][:, :])
close(fid)

function calculate_NICV(traces, variances, intermediate_value_vector)
    all_mean_values = CUDA.zeros(size(traces)[2], 256)
    for j in 0:255
        all_mean_values[:, j + 1] = mean(traces[intermediate_value_vector .== j, :], dims=1)[1, :]
    end
    return var(all_mean_values, dims=2) ./ variances
end

sections_of_trace = CuArray(original_matrix_of_current_data)
all_intermediate_values = CuArray(intermediate_value_vector)

variances = var(sections_of_trace, dims=1)[1, :]

nicv = calculate_NICV(sections_of_trace, variances, all_intermediate_values)

p = plot(Array(nicv))