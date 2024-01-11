function rotmat = rotx(alpha)
%rotx     Rotation matrix around x-axis
%   ROTMAT = rotx(ALPHA) returns the rotation matrix, ROTMAT, that rotates
%   a point around the x-axis for an angle ALPHA (in degrees). The point is
%   specified in the form of [x;y;z], with the x, y, and z axes forming a
%   right-handed Cartesian coordinate system. With the x axis pointing
%   towards the observer, ALPHA is measured counter-clockwise in the y-z
%   plane.
%
%   ROTMAT is a 3x3 matrix. The rotation of the point can be achieved by
%   left-multiplying ROTMAT with the point's coordinate vector [x;y;z].
%
%   % Example:
%   %   Rotate a point, (0,1,0), around x-axis 45 degrees
%   %   counter-clockwise.
%
%   p = [0;1;0];
%   p = rotx(45)*p

%   Copyright 2012 The MathWorks, Inc.

%   References:
%   [1] James Foley, et. al. Computer Graphics Principles and Practices in
%       C, 2nd Edition, Addison-Wesley, 1995

%#codegen
%#ok<*EMCA>

eml_assert_no_varsize(1,alpha);
sigdatatypes.validateAngle(alpha,'rotx','ALPHA',{'scalar'});
% rotate in the direction of y->z, counter-clockwise
rotmat = [1 0 0;0 cosd(alpha) -sind(alpha); 0 sind(alpha) cosd(alpha)];

% [EOF]
