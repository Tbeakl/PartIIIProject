using HDF5
include("../encryption/leakage_functions.jl")
# The parts of the trace we are interested in I belive are the outputs of the adds and xors (because the rotation will be heavily
# related to the output of the xor but for ease of putting into the factor graph I will put the values in after the rotate so we want to select with
# 10001)
elements_of_trace_to_select = append!(repeat([true, false, false, false, true], 320), ones(Bool, 16))



for file_number in 0:9
    intermediates_fid = h5open(string("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\intermediate_value_traces_2\\recording_attack_counter_constant_", file_number, ".hdf5"), "w")
    trace_fid = h5open(string("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\captures\\ChaChaRecordings_2\\recording_attack_counter_constant_", file_number, ".hdf5"), "r")

    for trace_number in 0:99
        for element_to_select in 0:14
            key = UInt32.(read(trace_fid[string("power_", trace_number, "_", element_to_select)]["key"]))
            nonce = UInt32.(read(trace_fid[string("power_", trace_number, "_", element_to_select)]["nonce"]))
            counter = UInt32.(read(trace_fid[string("power_", trace_number, "_", element_to_select)]["counter"]))[1]
            plaintext = UInt32.(read(trace_fid[string("power_", trace_number, "_", element_to_select)]["plaintext"]))
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
            intermediates_fid[string("power_", trace_number, "_", element_to_select)] = all_values
        end
    end

    close(trace_fid)
    close(intermediates_fid)
end
# Need to make the code for correctly selecting the parts of the trace we care about, and add the plaintext values into the output alongside
# the output of the ciphertext, key, nonce, counter



# Basic structure of the stored traces will be 
# |key|nonce|counter|plaintext|values_of_intermediates|ciphertext|