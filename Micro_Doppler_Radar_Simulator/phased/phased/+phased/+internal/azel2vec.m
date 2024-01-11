function [az_vec, el_vec] = azel2vec(azel)
%This function is for internal use only. It may be removed in the future.

%azel2vec   Compute the direction vector for azimuth and elevation
%   [AZ_VEC, EL_VEC] = phased.internal.azelvec(AZEL) computes the
%   increasing direction vectors, AZ_VEC and EL_VEC, at azimuth and
%   elevation angles (in degrees) specified in AZEL.
%
%   AZEL is a 2xN matrix whose columns are azimuth and elevation angle
%   pairs in the form of [az;el]. Azimuth angles must be between -180 and
%   180 while elevation angles must be between -90 and 90. AZ_VEC and
%   EL_VEC are 3xN matrix whose columns represent the corresponding
%   direction vectors in azimuth and elevation directions, respectively.
%   The vectors are in the form of [x;y;z].
%
%   % Example:
%   %   Determine the direction vectors at 0 degree azimuth and 45 degree
%   %   elevation, as well as 45 degrees azimuth and 0 degrees elevation.
%
%   [azvec, elvec] = phased.internal.azel2vec([0 45;45 0])

%   Copyright 2012 The MathWorks, Inc.

%#codegen

theta = azel(2,:);
phi = azel(1,:);
el_vec = [-sind(theta).*cosd(phi);...
    -sind(theta).*sind(phi);cosd(theta)];
az_vec = [-sind(phi);cosd(phi);zeros(1,numel(phi))];


% [EOF]
