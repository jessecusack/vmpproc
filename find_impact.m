function [idx] = find_impact(t_fast, P_fast, Ax, Ay, lat, Azc, Astdc, npoints)
% Parameters
% t_fast: fast time [s]
% P_fast: fast pressure [dbar]
% Ax: x acceleration [-]
% Ay: y acceleration [-]
% lat: latitude of profile [degrees north]
% Azc: Vertical acceleration threshold, default is 0.1 m s-2
% Astdc: Horizontal acceleration standard deviation threshold, default is
% 200 (in whatever units Rockland has devised for this varaible)
% npoints: Number of points for the moving standard deviation, default is 21
%
% Returns
% idx: index of last good data before bottom impact
%
% Notes
% * Requires the gsw toolbox. 
% * Assumes a downward profile where time values increase. 
% * Vertical acceleration is calculated as the second derivative of height with respect to time. 
% * Only looks at the last 1000 points of the profile, hopefully the data are not shorter than this!

addpath(genpath('../matlab_toolboxes/'))

if ~exist('Azc', 'var')
    Azc = 0.1;
end
if ~exist('Astdc', 'var')
    Astdc = 200;
end
if ~exist('npoints', 'var')
    npoints = 21;
end

% Only need the end of the profile, lets use the last 1000 points.
i1 = length(t_fast) - 1000;

z = gsw_z_from_p(P_fast(i1:end), lat);
Az = gradient(gradient(z, t_fast(i1:end)), t_fast(i1:end));

[~, imax] = max(abs(Az));

idx1 = find(abs(Az(1:imax)) < Azc, 1, 'last');

Ax_std = sqrt(movvar(Ax(i1:end), npoints));
Ay_std = sqrt(movvar(Ay(i1:end), npoints));
Axy_std_mean = 0.5*(Ax_std + Ay_std);

idx = i1 + find(abs(Axy_std_mean(1:idx1)) < Astdc, 1, 'last');

end