function M = azelcoord(az,el)
%This function is for internal use only. It may be removed in the future.

%azelcoord Coordinate system at azimuth and elevation direction
%   M = phased.internal.azelcoord(AZ,EL) returns the axes of coordinate
%   system, (az_hat, el_hat, r_hat), at the direction of (AZ,EL) (in
%   degrees) in the (x,y,z) system where az and el are defined. M is a 3x3
%   matrix whose columns represent the three axes, az_hat, el_hat, and
%   r_hat. AZ must be between -180 and 180. El must be between -90 and 90.
%
%   The coordinate system is constructed at the direction of (AZ,EL). The
%   three axes are az_hat, the direction of increasing azimuth; el_hat, the
%   direction of increasing elevation; and r_hat, the increasing radial
%   direction.
%
%   % Example:
%   %   Determine the local coordinate system defined at the direction of
%   %   45 degrees azimuth and 0 degrees elevation.
%   
%   M = phased.internal.azelcoord(45,0)

%   Copyright 2012 The MathWorks, Inc.

%#codegen

% az: x->y
% el: xy->z

% M = [...
%     -sin(az_rad) cos(az_rad) 0;...
%     -sin(el_rad)*cos(az_rad) -sin(el_rad)*sin(az_rad) cos(el_rad);...
%     cos(el_rad)*cos(az_rad) cos(el_rad)*sin(az_rad) sin(el_rad)].';

r_vec = [cosd(el).*cosd(az) cosd(el).*sind(az) sind(el)].';
[az_vec, el_vec] = phased.internal.azel2vec([az;el]);

M = [az_vec el_vec r_vec];



% [EOF]
