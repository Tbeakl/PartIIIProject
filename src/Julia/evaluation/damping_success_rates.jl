using HDF5

damping_factors = ["1_0", "0_99", "0_95", "0_9", "0_75", "0_5"]

for damping_factor in damping_factors
    fid = h5open(string("./data/evaluation/damping_experimentation/", damping_factor, "_results.hdf5"), "r")
    println(damping_factor, " ", mean(read(fid["successes"])))
    close(fid)
end