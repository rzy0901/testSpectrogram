function rotaxes = rotazel(startaxes,azel)
%This function is for internal use only. It may be removed in the future.

%ROTAZEL  Rotate axes for a given azimuth and elevation pair
%   AXES_NEW = phased.internal.rotazel(AXES_OLD,AZEL) rotates the axes
%   specified in AXES_OLD by a pair of azimuth/elevation angles specified
%   in AZEL (in degrees). AXES_OLD is a 3-row matrix with each column
%   specifying an axis in [x;y;z] form. Note that axis can be interpreted
%   as either a position vector (in meters) or a direction vector.
%
%   AZEL is a length-2 vector in the form of [azimuth; elevation] where
%   azimuth is within [-180 180] and elevation is within [-90 90]. AXES_NEW
%   has the same dimension as AXES_OLD, with each column representing the
%   rotated axis.
%
%   Example:
%   %   Rotate x axis by 90 degrees azimuth.
%   
%   axis_new = phased.internal.rotazel([1;0;0],[90;0])

%   Copyright 2011 The MathWorks, Inc.

%   Reference:
%   [1] http://en.wikipedia.org/wiki/Euler_angles, "Relationship with
%   Physical Motions"

% to make the rotation easier, we should first rotate with y and then
% rotate with z because rotating y does not change the rotation axis for
% rotating z. If we rotate z first, then we need to rotate with new y axis.

%#codegen
    
    rotaxes = rotationMatrixForZ(azel(1))*...
    rotationMatrixForY(azel(2))*startaxes;

end

function rotmat = rotationMatrixForZ(gamma)

% rotate in the direction of x->y
u = cosd(gamma);
v = sind(gamma);
rotmat = [u -v 0; v u 0; 0 0 1];

end

function rotmat = rotationMatrixForY(beta)

% rotate in the direction of x->z
u = cosd(beta);
w = sind(beta);
rotmat = [u 0 -w; 0 1 0; w 0 u];

end

% [EOF]
