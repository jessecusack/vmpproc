bpFiles = dir("./proc/*/binned_*.mat");
cbFile = "proc/all_binned.mat";

join_binned(bpFiles, cbFile, overwrite=true)
binned_to_netcdf(cbFile, overwrite=true)