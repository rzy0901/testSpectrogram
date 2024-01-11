function uv = azel2uv(azel)
%AZEL2UV  Convert angles from az/el format to u/v format
%   UV = azel2uv(AzEl) converts the azimuth/elevation angle pairs to their
%   corresponding u/v values. AzEl is a 2-row matrix whose columns are
%   angles specified in [azimuth;elevation] (in degrees) form. UV has the
%   same dimension as AzEl. The columns of UV are angles specified in [u;v]
%   format.
%
%   Assume the boresight is the X-axis. The azimuth angle is defined as the
%   angle from the X-axis toward the Y-axis, ranging from -90 to 90
%   degrees. The elevation angle is the angle from the X-Y plane toward the
%   Z-axis, ranging from -90 to 90 degrees. The phi angle is the angle from
%   the Y-axis to the Z-axis, ranging from 0 to 360 degrees. The theta
%   angle is the angle from the X-axis toward the Y-Z plane, ranging from 0
%   to 90 degrees. The value u is defined as sin(theta)*cos(phi) and the
%   value v is defined as sin(theta)*sin(phi).
%
%   % Example: 
%   %   Find the corresponding u and v representation for 30 degrees
%   %   azimuth and 0 degrees elevation.
%   
%   uv = azel2uv([30;0])
%
%   See also azel2phitheta, azel2uvpat, uv2azel, uv2azelpat.

%   Copyright 2011 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

sigdatatypes.validateAzElAngle(azel,'azel2uv','AzEl');

az = azel(1,:);
el = azel(2,:);

% az should be limited between [-90 90]
validateattributes(az,{'double'},{'>=',-90,'<=',90},...
    'azel2uv', 'angle az');

u = cosd(el).*sind(az);
v = sind(el);

uv = [u;v];

% [EOF]
