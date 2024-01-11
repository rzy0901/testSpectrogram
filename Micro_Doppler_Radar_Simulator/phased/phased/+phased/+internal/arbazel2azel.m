function azel = arbazel2azel(azel_in)
%This function is for internal use only. It may be removed in the future.

%ARBAZEL2AZEL   Convert arbitrary az/el angles to standard az/el angles
%   AZEL = phased.internal.arbazel2azel(AZEL_IN) converts the arbitrary
%   [azimuth; elevation] angles (in degrees) specified in AZEL_IN to
%   standard [azimuth; elevation] angles (in degrees) AZEL where azimuth is
%   within [-180 180] and elevation is within [-90 90]. AZEL_In is a 2-row
%   matrix whose columns represents az/el angles and AZEL has the same
%   dimension as AZEL_IN. Columns in AZEL contain the standard az/el angles
%   corresponding to the az/el angles in AZEL_IN.
%
%   % Example:
%   %   Convert 360 degrees azimuth and 120 degrees elevation into standard
%   %   azimuth elevation values.
%
%   azel = phased.internal.arbazel2azel([360;120])
%
%   See also phased, azel2dirvec, dirvec2azel

%   Copyright 2015 The MathWorks, Inc.

%#ok<*EMCA>
%#codegen

azel_in = phased.internal.deg2rad(azel_in);
[th,phi] = cart2sph(...
    cos(azel_in(1,:)).*cos(azel_in(2,:)),...
    sin(azel_in(1,:)).*cos(azel_in(2,:)),...
    sin(azel_in(2,:)));
azel = phased.internal.rad2deg([th;phi]);
