function rotmat = rotz(gamma)
%rotz     Rotation matrix around z-axis
%   ROTMAT = rotz(GAMMA) returns the rotation matrix, ROTMAT, that rotates
%   a point around the z-axis for an angle GAMMA (in degrees). The point is
%   specified in the form of [x;y;z], with the x, y, and z axes forming a
%   right-handed Cartesian coordinate system. With the z axis pointing
%   towards the observer, GAMMA is measured counter-clockwise in the x-y
%   plane.
%
%   ROTMAT is a 3x3 matrix. The rotation of the point can be achieved by
%   left-multiplying ROTMAT with the point's coordinate vector [x;y;z].
%
%   % Example:
%   %   Rotate a point, (0,1,0), around z-axis 45 degrees
%   %   counter-clockwise.
%
%   p = [0;1;0];
%   p = rotz(45)*p

%   Copyright 2012 The MathWorks, Inc.

%   References:
%   [1] James Foley, et. al. Computer Graphics Principles and Practices in
%       C, 2nd Edition, Addison-Wesley, 1995

%#codegen
%#ok<*EMCA>

eml_assert_no_varsize(1,gamma);
sigdatatypes.validateAngle(gamma,'rotz','GAMMA',{'scalar'});
% rotate in the direction of x->y, counter-clockwise
rotmat = [cosd(gamma) -sind(gamma) 0; sind(gamma) cosd(gamma) 0; 0 0 1];


% [EOF]
