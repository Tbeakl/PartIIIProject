using HDF5
include("../encryption/leakage_functions.jl")
# The parts of the trace we are interested in I belive are the outputs of the adds and xors (because the rotation will be heavily
# related to the output of the xor but for ease of putting into the factor graph I will put the values in after the rotate so we want to select with
# 10001)
elements_of_trace_to_select = append!(repeat([true, false, false, false, true], 320), ones(Bool, 16))

file_number = 1
trace_number = 1
# for file_number in 1:360
    # intermediates_fid = h5open(string("D:\\Year_4_Part_3\\Dissertation\\SourceCode\\PartIIIProject\\data\\intermediate_value_traces\\recording_profiling_", file_number, ".hdf5"), "w")
    trace_fid = h5open(string("D:\\ChaChaData\\captures\\ChaChaRecordings\\recording_profiling_", file_number, ".hdf5"), "r")


    # for trace_number in 0:249
        key = UInt32.(read(trace_fid[string("power_", trace_number)]["key"]))
        nonce = UInt32.(read(trace_fid[string("power_", trace_number)]["nonce"]))
        counter = UInt32.(read(trace_fid[string("power_", trace_number)]["counter"]))[1]
        plaintext = UInt32.(read(trace_fid[string("power_", trace_number)]["plaintext"]))
        ciphertext = encrypt(key, nonce, counter) .⊻ plaintext
        trace = encrypt_collect_trace(key, nonce, counter, byte_values_for_input)#[elements_of_trace_to_select]

        all_trace_values = append!(
            byte_values_for_input.(key),
            byte_values_for_input.(nonce),
            byte_values_for_input.([counter]),
            byte_values_for_input.(plaintext),
            trace,
            byte_values_for_input.(ciphertext))
        all_values = UInt8.(collect(Iterators.flatten(all_trace_values)))
        # intermediates_fid[string("power_", trace_number)] = all_values

        UInt8.(collect(Iterators.flatten(trace)))
    # end

    close(trace_fid)
    close(intermediates_fid)
# end
# Need to make the code for correctly selecting the parts of the trace we care about, and add the plaintext values into the output alongside
# the output of the ciphertext, key, nonce, counter



# Basic structure of the stored traces will be 
# |key|nonce|counter|plaintext|values_of_intermediates|ciphertext|