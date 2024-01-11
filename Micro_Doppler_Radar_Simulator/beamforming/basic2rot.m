% /*************************************************************************************
%    Intel Corp.
%
%    Project Name:  Conference Room Channel Model
%    File Name:     basic2rot.m
%    Authors:       A. Lomayev, R. Maslennikov
%    Version:       1.0
%    History:       May 2010 created
%
%  *************************************************************************************
%    Description:
%
%    function performs recalculation of azimuth and elevation angles
%    determinig in basic coordinates to rotated coordinates,
%    univocal rotated system position is set by 3 Euler's rotations
%    (by corresponding 3 angles)
%
%    [azr,elr] = basic2rot(az,el,az_rot,el_rot,self_rot)
%
%    Inputs:
%
%       1. az     - array of input azimuth angles in [deg] in basic coordinates
%       2. el     - array of input elevation angles in [deg] in basic coordinates
%       3. az_rot - angle in [deg] determines azimuth rotation for rotated system of coordinates relative to basic system
%       4. el_rot - angle in [deg] determines elevation rotation for rotated system of coordinates relative to basic system
%       5. sr_rot - angle in [deg] determines self rotation for rotated system of coordinates relative to basic system
%
%    Outputs:
%
%       1. az_rot - array of output azimuth angles in [deg] in rotated coordinates
%       2. el_rot - array of output elevation angles in [deg] in rotated coordinates
%
%  *************************************************************************************/
function [azr,elr] = basic2rot(az,el,az_rot,el_rot,self_rot)

% radian to degree and vice versa conversion constants
d2r = pi./180;
r2d = 180./pi;

% degree to radian conversion
az = az.*d2r;
el = el.*d2r;
az_rot = az_rot.*d2r;
el_rot = el_rot.*d2r;
self_rot = self_rot.*d2r;

% spherical to cartesian
x = sin(el).*cos(az);
y = sin(el).*sin(az);
z = cos(el);

% rotation matrix: 3 Euler's rotations
R_self = [cos(self_rot),sin(self_rot),0;-sin(self_rot),cos(self_rot),0;0,0,1];
R_el   = [1,0,0;0,cos(el_rot),sin(el_rot);0,-sin(el_rot),cos(el_rot)];
R_az   = [cos(az_rot),sin(az_rot),0;-sin(az_rot),cos(az_rot),0;0,0,1];
R = R_self*R_el*R_az;

% calculate coordinates in rotated system
dec_rot = R*[x.';y.';z.'];

xr = dec_rot(1,:);
yr = dec_rot(2,:);
zr = dec_rot(3,:);

% cartesian to spherical
azr = atan2(yr,xr); % range [-pi/2:pi/2]
elr = acos(zr);

% radian to degree
azr = azr.*r2d;
% azimuth angle should be [0:360]
azr = mod(360 + azr,360);
elr = elr.*r2d;