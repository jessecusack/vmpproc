addpath("./vmpproc/")
addpath("./odas/");
addpath(genpath("./gsw/"));

pFiles = dir("./raw/*.p");
saveDir = "./proc/";

% List of filenames to exclude from conductivity offset (these were
% calibration casts)
exclude = ["A142_0011.p", "A412_0026.p"];

% Processing parameters
info = get_info();
info.overwrite = false;
info.pMin = 1;

% Setup GPS file.
gps_file = "~/path/to/ship/data/ship.merged.nc";
time = ncread(gps_file, 't');
lon = ncread(gps_file, 'lon');
lat = ncread(gps_file, 'lat');

gps.time = time/86400 + datenum("1970-01-01");
gps.lon = lon;
gps.lat = lat;

% Extract profiles from p files 
for idx = 1:length(pFiles)
    name = pFiles(idx).name;

    if ismember(name, exclude)
        continue
    end

    pFile = fullfile(pFiles(idx).folder, name);

    % Account for conductivity offset in instrument 142
    if contains(name, 'A142')
        fprintf("Offsetting conductivity!")
        c_offset = 0.194;
    else
        c_offset = 0.0;
    end

    try
        generate_diss_profiles(pFile, saveDir, info=info, gps=gps, c_offset=c_offset);
    catch e
        fprintf(1, 'There was an error! The message was:\n%s\n', e.message);
    end
end