function azel = phitheta2azel(phitheta)
%PHITHETA2AZEL Convert angles from phi/theta format to az/el format
%   AzEl = phitheta2azel(PhiTheta) converts the phi/theta angle pairs to
%   their corresponding azimuth/elevation angles. PhiTheta is a 2-row
%   matrix whose columns are angles specified in [phi;theta] form. AzEl has
%   the same dimension as PhiTheta. The columns of AzEl are angles
%   specified in [azimuth;elevation] format. All angles are specified in
%   degrees.
%
%   Assume the boresight is the X-axis. The azimuth angle is defined as the
%   angle from the X-axis toward the Y-axis, ranging from -180 to 180
%   degrees. The elevation angle is the angle from the X-Y plane toward the
%   Z-axis, ranging from -90 to 90 degrees. The phi angle is the angle from
%   the Y-axis to the Z-axis, ranging from 0 and 360 degrees. The theta
%   angle is the angle from X-axis toward the Y-Z plane, ranging from 0 to
%   180 degrees.
%
%   % Example: 
%   %   Find the corresponding azimuth and elevation representation for 0 
%   %   degrees phi and 30 degrees theta.
%   
%   azel = phitheta2azel([0;30])
%
%   See also azel2phitheta, azel2phithetapat, phitheta2azelpat,
%   phitheta2uv.

%   Copyright 2011-2014 The MathWorks, Inc.

%#codegen 
%#ok<*EMCA>

validateattributes(phitheta,{'double'},{'2d','real','nonnegative'},...
    'phitheta2azel','PhiTheta');

cond = size(phitheta,1) == 2;
if ~cond
    coder.internal.assert(cond,'phased:phased:invalidRowNumbers','PhiTheta',2);
end
phi = phitheta(1,:);
theta = phitheta(2,:);

sigdatatypes.validateAngle(phi,'phitheta2azel','angle phi',...
    {'>=',0,'<=',360});
sigdatatypes.validateAngle(theta,'phitheta2azel','angle theta',...
    {'>=',0,'<=',180});


% azimuth 
az = atand(tand(theta).*cosd(phi));

% take care of quadrants
temp_idx = find((theta>90)&((phi<90)|(phi>270)));
az(temp_idx) = az(temp_idx)+180;
temp_idx = find((theta>90)&((phi>90)&(phi<270)));
az(temp_idx) = az(temp_idx)-180;

% take care of planes
az(((phi==90)|(phi==270))&(theta>90)) = 180;

% take care of axis
az(theta==180) = 180;
az((theta==90)&((phi==90)|(phi==270))) = 0;

% elevation
el = asind(sind(theta).*sind(phi));

azel = [az;el];




% [EOF]
