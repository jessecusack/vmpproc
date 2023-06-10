function [] = generate_diss_profiles(pFile, saveDir, varargin)
% This function consumes a p file and produces a individual profile mat
% files of dissipation and other quantities.
%
% Arguments
% ---------
% pFile : text
%    Data file to read. Must be RSI .p file. If a .mat file of the same
%    name exists with the .p file, it will be read preferentially. 
% saveDir : text
%    Directory to save the profiles. This function creates a subdirectory
%    with a name that matches the p file name.
% info : struct, optional (but highly recommended) parameter
%    Structure containing variables relevant to the dissipation
%    calculation. Similar to the ql_info structure returned by 
%    quick_look. See help get_info for the description of the structure and
%    see help for despike, get_profile or get_diss_odas to understand most
%    of the parameters. 
%
%    The fields in info that are unique to this function are:
%
%        chop_end : ["bottom", "top"], optional
%            Choose which end of the profile to chop data from in order to
%            fit an integer number of fft lengths. Default is "bottom". 
%            Choosing "top" is optimal for estimating bottom boundary layer
%            dissipation.
%        overwrite : [true, false], optional
%            Choose whether to overwrite existing profiles. Default is
%            false.
%        revert_to_defalt_gps [true, false], optional
%            If false, an error will be thrown if the GPS time series does 
%            not cover all the profile times. If true, the error is reduced 
%            to a warning. Default is false.
%
% gps : struct, optional 
%     Structure containing the `lon` and `lat`. Optionally is may contain 
%     the field time. Time must be in matlab datenum type and the same 
%     length as `lon` and `lat`. The fields may have length 1, or be a 
%     series. mIn the latter case, coordinates are linearly interpolated 
%     to profile time.
% c_offset : number, optional 
%     Offset applied to JAC conductivity prior to any calculations. A basic 
%     method of calibration. Default is 0.0 mS cm-1. 
%
% First created by Jesse Cusack (jesse.cusack@oregonstate.edu) 2023-06-02.

lon_default = -45;
lat_default = 45;

% Parse arguments
iP = inputParser;
iP.StructExpand = false;
validText = @(x) isstring(x) || ischar(x);
validNumber = @(x) isnumeric(x) && isscalar(x);
addRequired(iP,'pFile', validText);
addRequired(iP,'saveDir', validText);
addParameter(iP,'info', get_info(), @isstruct);
addParameter(iP, 'gps', struct('lon', lon_default, 'lat', lat_default, 'default', true), @isstruct)
addParameter(iP, 'c_offset', 0, validNumber)
parse(iP, pFile, saveDir, varargin{:}); 
pFile = iP.Results.pFile;
saveDir = iP.Results.saveDir;
info = iP.Results.info;
gps = iP.Results.gps;
c_offset = iP.Results.c_offset;

fprintf("\nP file: %s\n", pFile)

[pPath, name, ~] = fileparts(pFile);
% Check if the directory exists
savePath = fullfile(saveDir, name);
if exist(savePath, "dir")
    if (length(dir(fullfile(savePath, "profile_*.mat"))) > 1) && ~info.overwrite
        error("Save path (%s) exists and contains profiles. Overwrite is false.", savePath)
    end
end

matFile = fullfile(pPath, strcat(name, ".mat"));

p = odas_p2mat(pFile);

if ~exist(matFile, "file")
    fprintf("\nSaving p file to mat.\n")
    save(matFile, '-struct', 'p')
end

% Extract serial number from setup file
sn = str2double(regexp(p.setupfilestr, '(?<=sn\s*=[^0-9]*)[0-9]*\.?[0-9]+', 'match'));

% Find profiles
profileIdx = get_profile(p.P_fast, p.W_fast, info.pMin, info.wMin, ...
    'down', info.minDuration, p.fs_fast);
nProfiles = size(profileIdx, 2);

if nProfiles < 1
    error("No profiles detected.")
end

profileIdxSlow = zeros(size(profileIdx));
for i = 1:nProfiles
    i1 = profileIdx(1, i);
    i2 = profileIdx(2, i);
%     ii = i1 + find_impact(p.t_fast(i1:i2), p.P_fast(i1:i2), p.Ax(i1:i2), p.Ay(i1:i2), lat, 0.1, 200, 21);
    profileIdxSlow(1, i) = find(p.t_slow > p.t_fast(i1), 1, 'first');
    profileIdxSlow(2, i) = find(p.t_slow < p.t_fast(i2), 1, 'last');
end

time_start = p.filetime + p.t_slow(profileIdxSlow(1, :))/86400;
time_end = p.filetime + p.t_slow(profileIdxSlow(2, :))/86400;

% Check GPS
[default_gps, interp_gps] = check_gps(gps);
gps_overlap = true(nProfiles, 1);
if interp_gps
    gps = clean_gps(gps);
    gps_overlap = (gps.time(1) < time_start) & (gps.time(end) > time_start);
    msg = sprintf("GPS times overlaps with %i of %i profiles.", sum(gps_overlap), nProfiles);
    if ~all(gps_overlap) && info.revert_to_default_gps
        warning(msg)
    elseif ~all(gps_overlap) && ~info.revert_to_default_gps
        error(msg)
    else
        fprintf(strcat(msg, "\n"))
    end
end

fprintf("\nCalculating dissipation with the following parameters:\n")
disp(info)

% Create save directory
mkdir(saveDir, name)
fprintf("\nSaving profiles to %s\n", savePath)

% Put in info variables
info.fs_fast = p.fs_fast;
info.fs_slow = p.fs_slow;

% Offset conductivity
p.JAC_C = p.JAC_C + c_offset;

for idx = 1:nProfiles
    fprintf('Profile %i/%i\n', idx, nProfiles)


    dt = datestr(time_start(idx), "yyyymmddTHHMMSSZ");
    saveName = sprintf("profile_%s.mat", dt);
    saveFullFile = fullfile(savePath, saveName);
%     if exist(saveFullFile, "file") && ~info.overwrite
%         fprintf("%s already exists, skipping.\n", saveName)
%         continue
%     end

    pfl = struct;
    pfl.sn = sn;
    
    idx_start = profileIdxSlow(1, idx);
    idx_end = profileIdxSlow(2, idx);

    idxs_slow = idx_start:idx_end;

%     % THIS HERE IS THE KEY BIT FOR BOTTOM DISSIPATION!!!
%     % restrict range to capture bottom
%     idx0 = mod(length(SH1_des), diss_info.diss_length) + 1;
%     range = idx0:length(SH1_des);
%     ninrange = length(range);

    idxs_fast = profileIdx(1, idx):profileIdx(2, idx);

    pfl.time_start = time_start(idx);
    pfl.time_end = time_end(idx);
    
    % attach GPS
    if interp_gps && gps_overlap(idx)
        lon = interp1(gps.time, gps.lon, pfl.time_start);
        lat = interp1(gps.time, gps.lat, pfl.time_start);
        pfl.lon = lon;
        pfl.lat = lat;
    elseif (interp_gps && ~gps_overlap(idx)) || default_gps
        fprintf("Using default GPS coordinates (Lon, Lat) = (%i, %i).\n", ...
            lon_default, lat_default)
        lon = lon_default;
        lat = lat_default;
        pfl.lon = NaN;
        pfl.lat = NaN;
    elseif ~interp_gps && ~default_gps
        lon = gps.lon;
        lat = gps.lat;
        pfl.lon = lon;
        pfl.lat = lat;
    end

    % DESPIKE
    % Piezo-accelerometers
    [pfl.Ax_ds, ~, ~, pfl.Ax_frac] = despike(p.Ax(idxs_fast), info.thresh, info.smooth, p.fs_fast, round(info.N*p.fs_fast));
    [pfl.Ay_ds, ~, ~, pfl.Ay_frac] = despike(p.Ay(idxs_fast), info.thresh, info.smooth, p.fs_fast, round(info.N*p.fs_fast));
    % Shear probes
    [pfl.sh1_ds, ~, ~, pfl.sh1_frac] = despike(p.sh1(idxs_fast), info.thresh, info.smooth , p.fs_fast, round(info.N*p.fs_fast));
    [pfl.sh2_ds, ~, ~, pfl.sh2_frac] = despike(p.sh2(idxs_fast), info.thresh, info.smooth , p.fs_fast, round(info.N*p.fs_fast));
    % Conductivity
    [pfl.C, ~, ~, ~] = despike(p.JAC_C(idxs_slow), info.thresh, info.smooth , p.fs_slow, round(info.N*p.fs_slow));
    
    % The high pass frequency is half the fft length
    fc = 0.5/info.fft_length;
    [bh, ah] = butter(1, fc/(p.fs_fast/2), 'high');
    pfl.sh1_hp = filtfilt(bh, ah, pfl.sh1_ds);
    pfl.sh2_hp = filtfilt(bh, ah, pfl.sh2_ds);

    % DISSIPATION CALCULATION
    info.T = interp1(p.t_slow(idxs_slow), p.JAC_T(idxs_slow), p.t_fast(idxs_fast));
    info.t = p.t_fast(idxs_fast);
    info.P = p.P_fast(idxs_fast);

    % Do my own speed calculation because I'm not sure what they are doing
    info.speed = abs(gradient(gsw_z_from_p(info.P, lat), info.t));
    pfl.speed = info.speed;
    
    AA = [pfl.Ax_ds, pfl.Ay_ds]; 
    diss = get_diss_odas([pfl.sh1_hp, pfl.sh2_hp], AA, info);

    % FILL STRUCTURE
    % Slow variables
    pfl.P_slow = p.P_slow(idxs_slow);
    pfl.z_slow = gsw_z_from_p(pfl.P_slow, lat);
    pfl.T = p.JAC_T(idxs_slow);

    process_info = struct();
    process_info.speed = mean(pfl.speed);
    process_info.fs = round(p.fs_slow);
    pfl.SP = salinity_JAC_gsw(pfl.P_slow, pfl.T, pfl.C, process_info);
%     pfl.JAC_S = salinity_JAC(pfl.P_slow, pfl.T, pfl.C, process_info);

    pfl.t_slow = p.t_slow(idxs_slow);
    pfl.dn_slow = pfl.t_slow/86400 + p.filetime;
    % Thermodynamics
    pfl.SA = gsw_SA_from_SP(pfl.SP, pfl.P_slow, lon, lat);
    pfl.CT = gsw_CT_from_t(pfl.SA, pfl.T, pfl.P_slow);
%     [N2, P_mid] = gsw_Nsquared(pfl.SA, pfl.CT, pfl.P_slow, lat);
%     pfl.N2 = interp1(P_mid, N2, pfl.P_slow);
    pfl.rho0 = gsw_pot_rho_t_exact(pfl.SA, pfl.T, pfl.P_slow, 0);

    % Epsilon
    pfl.P_diss = diss.P;
    pfl.z_diss = gsw_z_from_p(pfl.P_diss, lat);
    pfl.eps1 = diss.e(1, :)';
    pfl.eps2 = diss.e(2, :)';
    
    % fast variables
    pfl.P_fast = info.P;
    pfl.z_fast = gsw_z_from_p(pfl.P_fast, lat);
    pfl.t_fast = info.t;
    pfl.dn_fast = pfl.t_fast/86400 + p.filetime;
    pfl.T1_fast = p.T1_fast(idxs_fast);
    pfl.T2_fast = p.T2_fast(idxs_fast);

    % SAVE
    fprintf("Saving to %s\n", saveFullFile)
    save(saveFullFile, "-struct", "pfl")

end

end