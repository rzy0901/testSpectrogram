function dirvec = azel2dirvec(azel)
%This function is for internal use only. It may be removed in the future.

%AZEL2DIRVEC    Convert azimuth/elevation angles to directional vectors
%   DIRVEC = phased.internal.azel2dirvec(AZEL) converts the angle specified
%   in AZEL of form [azimuth; elevation] (in degrees) to the corresponding
%   directional vector DIRVEC. AZEL is a 2-row matrix and DIRVEC is a 3-row
%   matrix whose columns contain the directional vectors of corresponding
%   angles in AZEL.
%
%   The directional vector is defined as [cos(el)cos(az); cos(el)sin(az);
%   sin(el)].
%
%   % Example:
%   %   Compute the directional vector corresponding to 45 degrees azimuth 
%   %   and 45 degrees elevation.
%
%   dirvec = phased.internal.azel2dirvec([45;45])
%
%   See also phased, arbazel2azel, dirvec2azel

%   Copyright 2015 The MathWorks, Inc.

%#ok<*EMCA>
%#codegen

azel = phased.internal.deg2rad(azel);
dirvec = [cos(azel(2,:)).*cos(azel(1,:));...
    cos(azel(2,:)).*sin(azel(1,:));sin(azel(2,:))];

