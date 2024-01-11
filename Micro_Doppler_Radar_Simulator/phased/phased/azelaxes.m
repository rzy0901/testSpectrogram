function ax = azelaxes(az,el)
%AZELAXES   Axes at given azimuth and elevation direction
%   AX = azelaxes(AZ,EL) returns the components of the set of three
%   orthonormal basis vectors at a point on the unit sphere. The basis
%   vectors are unit vectors in the radial, azimuthal, and elevation
%   directions. The point on the sphere is specified by the azimuth angle
%   AZ (in degrees) and the elevation angle EL (in degrees).
%
%   AX is a 3x3 matrix whose columns contain the unit vectors in the
%   radial, azimuthal, and elevation directions, respectively.
%
%   % Example:
%   %   Determine the basis vectors at a point on the sphere given by 30 
%   %   degrees azimuth and 0 degrees elevation.
%
%   ax = azelaxes(30,0)

%   Copyright 2012 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

eml_assert_no_varsize(1:nargin, az,el);
sigdatatypes.validateAngle(az,'azelaxes','AZ',{'scalar','<=',180,'>=',-180});
sigdatatypes.validateAngle(el,'azelaxes','EL',{'scalar','<=',90,'>=',-90});

ax = phased.internal.azelcoord(az,el);
ax = ax(:,[3 1 2]);

% [EOF]
