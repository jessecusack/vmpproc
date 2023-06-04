function [gps] = clean_gps(gps)
% Remove non-finte values
finite = isfinite(gps.time) & isfinite(gps.lon) & isfinite(gps.lat);
gps.time = gps.time(finite);
gps.lon = gps.lon(finite);
gps.lat = gps.lat(finite);
% Make unique
[~, idx, ~] = unique(gps.time);
gps.time = gps.time(idx);
gps.lon = gps.lon(idx);
gps.lat = gps.lat(idx);
end