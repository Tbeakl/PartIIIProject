using Random, HDF5
Random.seed!(1234)
include("../encryption/leakage_functions.jl")

random_keys::Vector{Vector{UInt32}} = []
random_nonces::Vector{Vector{UInt32}} = []
random_counters::Vector{UInt32} = []

fid = h5open("./PartIIIProject/data/evaluation/key_sets/100_random_keys_nonces_counter.hdf5", "w")
for i in 1:100
    fid[string("key_", i)] = generate_random_key()
    fid[string("nonce_", i)] = generate_random_nonce()
    fid[string("counter_", i)] = generate_random_counter()
end
close(fid)