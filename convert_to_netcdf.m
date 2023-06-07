function convert_to_netcdf(binnedFile, varargin)
% Convert a binned filed into a netcdf file. The netcdf is saved with the
% same filename as the mat file.
%
% Arguments
% ---------
% binnedFile : text
%    Data file to read.
% overwrite : [true false], optional
%    If true, overwrite existing netcdf file. Default is false.

% Parse arguments
iP = inputParser;
validText = @(x) isstring(x) || ischar(x);
addRequired(iP,'binnedFile', validText);
addParameter(iP,'overwrite', false, @islogical);
parse(iP, binnedFile, varargin{:}); 
overwrite = iP.Results.overwrite;

[path, name, ~] = fileparts(binnedFile);
ncFile = fullfile(path, strcat(name, ".nc"));

if exist(ncFile, "file") && ~overwrite
    error("%s exists, stopping.", ncFile)
else
    bd = load(binnedFile);
    delete(ncFile)
end

nz = length(bd.z);
nt = length(bd.time_start);

% Create coordinates
nccreate(ncFile, "z", "Dimensions", {"z" nz}, "Datatype", "single");
nccreate(ncFile, "time", "Dimensions", {"time" nt})
% Create data variables
nccreate(ncFile, "sn", "Dimensions", {"time" nt}, "Datatype", "int16")
nccreate(ncFile, "lon", "Dimensions", {"time" nt})
nccreate(ncFile, "lat", "Dimensions", {"time" nt})
nccreate(ncFile, "cast", "Dimensions", {"time" nt}, "Datatype", "int16")
nccreate(ncFile, "T", "Dimensions", {"z" nz "time" nt})
nccreate(ncFile, "SP", "Dimensions", {"z" nz "time" nt})
nccreate(ncFile, "CT", "Dimensions", {"z" nz "time" nt})
nccreate(ncFile, "SA", "Dimensions", {"z" nz "time" nt})
nccreate(ncFile, "rho0", "Dimensions", {"z" nz "time" nt})
nccreate(ncFile, "eps1", "Dimensions", {"z" nz "time" nt})
nccreate(ncFile, "eps2", "Dimensions", {"z" nz "time" nt})

% Write coordinates
ncwrite(ncFile, "z", bd.z);
ncwrite(ncFile, "time", (bd.time_start - datenum("1970-01-01"))*86400);  % Convert to posix time
% Write data variables
ncwrite(ncFile, "sn", bd.sn);
ncwrite(ncFile, "cast", bd.cast);
ncwrite(ncFile, "lon", bd.lon);
ncwrite(ncFile, "lat", bd.lat);
ncwrite(ncFile, "T", bd.T);
ncwrite(ncFile, "SP", bd.SP);
ncwrite(ncFile, "CT", bd.CT);
ncwrite(ncFile, "SA", bd.SA);
ncwrite(ncFile, "rho0", bd.rho0);
ncwrite(ncFile, "eps1", bd.eps1);
ncwrite(ncFile, "eps2", bd.eps2);

% Write attributes
addatt(ncFile, "z", "height", "Height", "m")
ncwriteatt(ncFile, "z", "positive", "up")
addatt(ncFile, "time", "time", "Time", "seconds since 1970-01-01")
ncwriteatt(ncFile, "sn", "long_name", "Serial number")
ncwriteatt(ncFile, "cast", "long_name", "Cast number")
addatt(ncFile, "lon", "longitude", "Longitude", "degree_east")
addatt(ncFile, "lat", "latitude", "Latitude", "degree_north")
addatt(ncFile, "T", "sea_water_temperature", "Temperature", "degree_C")
addatt(ncFile, "SP", "sea_water_practical_salinity", "Practical salinity", "")
addatt(ncFile, "CT", "sea_water_conservative_temperature", "Conservative temperature", "degree_C")
addatt(ncFile, "SA", "sea_water_absolute_salinity", "Absolute salinity", "g kg-1")
addatt(ncFile, "rho0", "sea_water_potential_density", "Potential density", "kg m-3")
ncwriteatt(ncFile, "rho0", "reference_pressure", "0 dbar")
addatt(ncFile, "eps1", "specific_turbulent_kinetic_energy_dissipation_in_sea_water", "Epsilon", "W kg-1")
addatt(ncFile, "eps2", "specific_turbulent_kinetic_energy_dissipation_in_sea_water", "Epsilon", "W kg-1")
end


function addatt(fn, loc, sn, ln, u)
ncwriteatt(fn, loc, "standard_name", sn)
ncwriteatt(fn, loc, "long_name", ln)
ncwriteatt(fn, loc, "units", u)
end


