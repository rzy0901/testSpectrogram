function azel = uv2azel(uv)
%UV2AZEL  Convert angles from u/v format to az/el format
%   AZEL = uv2azel(UV) converts the u/v value pairs to their corresponding
%   azimuth/elevation angles. UV is a 2-row matrix whose columns are values
%   specified in [u;v] format, where both u and v are between -1 and 1.
%   AzEl has the same dimension as UV. The columns of AzEl are the angles
%   specified in [azimuth;elevation] (in degrees) form.
%
%   Assume the boresight is the X-axis. The azimuth angle is defined as the
%   angle from the X-axis toward the Y-axis, ranging from -90 to 90
%   degrees. The elevation angle is the angle from the X-Y plane toward the
%   Z-axis, ranging from -90 to 90 degrees. The phi angle is the angle from
%   the Y-axis to the Z-axis, ranging from 0 and 360 degrees. The theta
%   angle is the angle from X-axis toward the Y-Z plane, ranging from 0 to
%   90 degrees. The value u is defined as sin(theta)*cos(phi) and the value
%   v is defined as sin(theta)*sin(phi).
%
%   % Example: 
%   %   Find the corresponding azimuth and elevation angles for u of 0.5 
%   %   and v of 0.
%   
%   azel = uv2azel([0.5;0])
%
%   See also azel2uv, azel2uvpat, uv2azelpat, uv2phitheta.

%   Copyright 2011-2014 The MathWorks, Inc.

%#codegen 
%#ok<*EMCA>


% u and v between [-1 1], within unit circle
validateattributes(uv,{'double'},{'2d','real','>=',-1,'<=',1},...
    'uv2azel','UV');

cond = size(uv,1) == 2;
if ~cond
    coder.internal.assert(cond,'phased:phased:invalidRowNumbers','UV',2);
end
u = uv(1,:);
v = uv(2,:);

cond = all(hypot(u,v)<=1);
if ~cond
    coder.internal.assert(cond,'phased:phased:invalidUV');
end

% az
temp = u.^2+v.^2;
temp(abs(temp)>1)=1;
az = atand(u./sqrt(1-temp));

% take care of quadrants

% take care of planes

% take care of axes
az((u==0)&(abs(v)==1))=0;

% el
el = asind(v);

azel = [az;el];

% [EOF]
