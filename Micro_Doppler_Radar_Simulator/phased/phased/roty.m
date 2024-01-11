function rotmat = roty(beta)
%roty     Rotation matrix around y-axis
%   ROTMAT = roty(BETA) returns the rotation matrix, ROTMAT, that rotates
%   a point around the y-axis for an angle BETA (in degrees). The point is
%   specified in the form of [x;y;z], with the x, y, and z axes forming a
%   right-handed Cartesian coordinate system. With the y axis pointing
%   towards the observer, BETA is measured counter-clockwise in the z-x
%   plane.
%
%   ROTMAT is a 3x3 matrix. The rotation of the point can be achieved by
%   left-multiplying ROTMAT with the point's coordinate vector [x;y;z].
%
%   % Example:
%   %   Rotate a point, (0,1,0), around y-axis 45 degrees
%   %   counter-clockwise.
%
%   p = [1;0;0];
%   p = roty(45)*p

%   Copyright 2012 The MathWorks, Inc.

%   References:
%   [1] James Foley, et. al. Computer Graphics Principles and Practices in
%       C, 2nd Edition, Addison-Wesley, 1995


%#codegen
%#ok<*EMCA>

eml_assert_no_varsize(1,beta);
sigdatatypes.validateAngle(beta,'roty','BETA',{'scalar'});
% rotate in the direction of z->x, counter-clockwise
rotmat = [cosd(beta) 0 sind(beta); 0 1 0; -sind(beta) 0 cosd(beta)];



% [EOF]
