function generate_binned_profiles(pflFiles, varargin)
% This bins a list of profile files. The minimum and maximum bin depths are
% found automatically. The binned profiles are saved into the same
% directories as the profile files. 
%
% Arguments
% ---------
% pFiles : struct
%     Structure generated by global search with dir. 
%     E.g. dir("folder/profile_*.mat")
% dz : numeric, optional
%     Bin spacing in m. Default is 2.
% zmin : numeric, optional
%     Minimum height bin (max depth) e.g. -200 m. Must be a negative 
%     number. If unspecified, zmin is determined automatically. 
% zmax : numeric, optional
%     Minimum height bin (min depth) e.g. -5 m. Must be a negative 
%     number. If unspecified, zmax is determined automatically.
% overwrite : bool, optional
%     Specify overwrite=true to overwrite existing binned files. Default is
%     false. 
%
% First created by Jesse Cusack (jesse.cusack@oregonstate.edu) 2023-06-04.

zmin_default = -999999;
zmax_default = -999999;

iP = inputParser;
iP.StructExpand = false;
validNumber = @(x) isnumeric(x) && isscalar(x);
validNegativeNumber = @(x) validNumber(x) && (x <= 0);
addRequired(iP, 'pflFiles', @isstruct);
addParameter(iP, 'dz', 2, validNumber);
addParameter(iP, 'zmin', zmin_default, validNegativeNumber);
addParameter(iP, 'zmax', zmax_default, validNegativeNumber);
addParameter(iP, 'overwrite', false);
parse(iP, pflFiles, varargin{:}); 
pflFiles = iP.Results.pflFiles;
dz = iP.Results.dz;
zmin = iP.Results.zmin;
zmax = iP.Results.zmax;
overwrite = iP.Results.overwrite;

nFiles = length(pflFiles);

% Find min/max height
find_zmin = zmin == zmin_default;
find_zmax = zmax == zmax_default;

if find_zmin && find_zmax
    zmin = 0;
    zmax = -100;
    for idx = 1:nFiles
        pflFile = fullfile(pflFiles(idx).folder, pflFiles(idx).name);
        pfl = load(pflFile, "z_slow");
        zmin = min([zmin; pfl.z_slow]);
        zmax = max([zmax; pfl.z_slow]);
    end
elseif find_zmin && ~find_zmax
    zmin = 0;
    for idx = 1:nFiles
        pflFile = fullfile(pflFiles(idx).folder, pflFiles(idx).name);
        pfl = load(pflFile, "z_slow");
        zmin = min([zmin; pfl.z_slow]);
    end
elseif ~find_zmin && find_zmax
    zmax = -100;
    for idx = 1:nFiles
        pflFile = fullfile(pflFiles(idx).folder, pflFiles(idx).name);
        pfl = load(pflFile, "z_slow");
        zmax = max([zmax; pfl.z_slow]);
    end
end

zmin = floor(zmin);
zmax = ceil(zmax);

for idx = 1:nFiles
    fprintf('Binning %i/%i\n', idx, nFiles)

    saveName = sprintf("binned_%s.mat", pflFiles(idx).name(9:11));
    saveFullFile = fullfile(pflFiles(idx).folder, saveName);
    if exist(saveFullFile, "file") && ~overwrite
        fprintf("%s already exists, skipping.\n", saveFullFile)
        continue
    end

    pflFile = fullfile(pflFiles(idx).folder, pflFiles(idx).name);
    pfl = load(pflFile);
    bp = bin_profile(pfl, zmin, zmax, dz);
    save(saveFullFile, "-struct", "bp")
end

end