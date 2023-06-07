# VMP processing toolbox

## Requirements

* `odas` - a proprietary matlab package written by Rockland Scientific Instruments
* `gsw` - the latest Gibbs Sea Water thermodynamic toolbox

## Usage

Processing follows three steps:
1. Convert Rockland's raw `.p` file format into a series of `.mat` profiles using `generate_diss_profiles.m`.
2. Depth bin average the generated profiles using `generate_binned_profiles.m`.
3. Combine the binned profiles into one dataset using `join_binned.m`.

Optionally, convert the joined dataset into a netCDF file using `convert_to_netcdf.m`.
