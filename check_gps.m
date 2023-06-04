function [default_gps, interp_gps] = check_gps(gps)

if isfield(gps, "default")
    fprintf("No gps file detected.\n" + ...
        "Using default coordinates (Lon, Lat) = (%i, %i).\n", ...
        gps.lon, gps.lat);
    default_gps = true;
    interp_gps = false;
elseif isfield(gps, "lon") && isfield(gps, "lat") && isfield(gps, "time")
    default_gps = false;
    llon = length(gps.lon);
    llat = length(gps.lat);
    lt = length(gps.time);
    if (llon == llat) && (llon == lt) && (lt > 1)
        interp_gps = true;
    elseif (llon == llat) && (llon == lt) && (lt == 1)
        interp_gps = false;
    else
        error("gps fields `time`, `lon`, and `lat` are not equal length.")
    end
elseif isfield(gps, "lon") && isfield(gps, "lat")
    if (length(gps.lon) == length(gps.lat)) && (length(gps.lon) == 1)
        interp_gps = false;
    else
        error("gps fields `lon` and `lat` are not both of length 1 or time is not specified.")
    end
else
    error("gps structure must contain fields `lon` and `lat`, and optionally `time`.") 
end

end