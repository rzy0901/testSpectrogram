function phitheta = azel2phitheta(azel)
%AZEL2PHITHETA Convert angles from az/el format to phi/theta format
%   PhiTheta = azel2phitheta(AzEl) converts the azimuth/elevation angle
%   pairs to their corresponding phi/theta angles. AzEl is a 2-row matrix
%   whose columns are angles specified in [azimuth;elevation] form.
%   PhiTheta has the same dimensions as AzEl. The columns of PhiTheta are
%   angles specified in [phi;theta] format. All angles are specified in
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
%   %   Find the corresponding phi and theta representation for 30 degrees
%   %   azimuth and 0 degrees elevation.
%   
%   phitheta = azel2phitheta([30;0])
%
%   See also azel2phithetapat, azel2uv, phitheta2azel, phitheta2azelpat.

%   Copyright 2011 The MathWorks, Inc.

%#codegen

sigdatatypes.validateAzElAngle(azel,'azel2phitheta','AzEl');

az = azel(1,:);
el = azel(2,:);

% phi
phi = atand(tand(el)./sind(az));

% take care of quadrants
temp_idx = find((az<0));
phi(temp_idx) = phi(temp_idx)+180;
temp_idx = find((az>0)&(el<0));
phi(temp_idx) = phi(temp_idx)+360;

% take care of planes
phi(((az==0)|(az==180)|(az==-180)) & (el>0)) = 90; 
phi(((az==0)|(az==180)|(az==-180)) & (el<0)) = 270; 

% take care of axes
phi(((az==0)|(az==180)|(az==-180)) & (el==0)) = 0; 
phi((az==90)&(el==0)) = 0;
phi((az==-90)&(el==0)) = 180;
phi(el==90) = 90;
phi(el==-90) = 270;

% theta
theta = acosd(cosd(el).*cosd(az));

phitheta = [phi;theta];

% [EOF]
