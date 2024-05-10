using HDF5
include("../encryption/leakage_functions.jl")
# The parts of the trace we are interested in I belive are the outputs of the adds and xors (because the rotation will be heavily
# related to the output of the xor but for ease of putting into the factor graph I will put the values in after the rotate so we want to select with
# 10001)
elements_of_trace_to_select = append!(repeat([true, false, false, false, true], 320), ones(Bool, 16))

path_to_data = "C:/Users/henry/Documents/PartIIIProject/data/"

# file_number = 1
# trace_number = 1
for file_number in 35:64
    println(file_number)
    intermediates_fid = h5open(string(path_to_data, "intermediate_value_traces_8_on_32/recording_profiling_", file_number, ".hdf5"), "w")
    trace_fid = h5open(string(path_to_data, "captures/ChaChaRecordings_8_on_32/recording_profiling_", file_number, ".hdf5"), "r")

    for trace_number in 0:999
        key = UInt32.(read(trace_fid[string("power_", trace_number)]["key"]))
        nonce = UInt32.(read(trace_fid[string("power_", trace_number)]["nonce"]))
        counter = UInt32.(read(trace_fid[string("power_", trace_number)]["counter"]))[1]
        plaintext = UInt32.(read(trace_fid[string("power_", trace_number)]["plaintext"]))
        ciphertext = encrypt(key, nonce, counter) .⊻ plaintext
        trace = encrypt_collect_trace(key, nonce, counter, full_value)[elements_of_trace_to_select]

        all_trace_values = append!(
            key,
            nonce,
            [counter],
            plaintext,
            trace,
            ciphertext)
        all_values = UInt32.(collect(Iterators.flatten(all_trace_values)))
        intermediates_fid[string("power_", trace_number)] = all_values
    end

    close(trace_fid)
    close(intermediates_fid)
end
# Need to make the code for correctly selecting the parts of the trace we care about, and add the plaintext values into the output alongside
# the output of the ciphertext, key, nonce, counter



# Basic structure of the stored traces will be 
# |key|nonce|counter|plaintext|values_of_intermediates|ciphertext|