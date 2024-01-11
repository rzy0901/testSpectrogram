function Ar = sph2cartvec(As,az,el)
%sph2cartvec Convert vector from spherical to Cartesian representation
%   VR = sph2cartvec(VS,AZ,EL) converts the vector VS from its spherical
%   representation to its corresponding Cartesian representation, VR. 
%
%   The spherical representation of a vector is its projection into a local
%   orthonormal basis defined by increasing values of azimuth angle
%   (az_hat), elevation angle (el_hat), and radius (r_hat). The basis
%   varies according its location on a sphere which is determined by the
%   azimuth and elevation angles, AZ and EL (in degrees). AZ must be
%   between -180 and 180 and EL must be between -90 and 90.
%
%   VS is a 3-row matrix whose columns are vectors in the spherical
%   representation, in the form of [az_hat;el_hat;r_hat]. VR has the same
%   dimensions as VS. Each column of VR contains the Cartesian
%   representation, in the form of [x;y;z], of the corresponding column of
%   VS.
%
%   % Example:
%   %   Determine the Cartesian representation of a vector, [1;1;1],
%   %   defined in the spherical coordinate at 30 degrees azimuth and 0
%   %   degrees elevation.
%
%   vr = sph2cartvec([1;1;1],30,0)

%   Copyright 2012 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

%   translates a vector in coordinates (az_hat, el_hat, r_hat), defined
%   in the spherical system at direction (az;el) to its corresponding
%   rectangular coordinates.

% Ar is 3xN rectangular coordinates [x; y; z]
% As is 3xN spherical coordinates [az;el;r] (see cart2sph)
% az is x->y in degrees, scalar
% el is xy->z in degrees, scalar

eml_assert_no_varsize(2:3,As,az,el);
validateattributes(As,{'double'},{'finite','nonnan','nonempty',...
    '2d','nrows',3},'sph2cartvec','VS');
sigdatatypes.validateAngle(az,'sph2cartvec','AZ',...
    {'scalar','<=',180,'>=',-180});
sigdatatypes.validateAngle(el,'sph2cartvec','EL',...
    {'scalar','<=',90,'>=',-90});

M = phased.internal.azelcoord(az,el);
    
Ar = M*As;


% [EOF]
