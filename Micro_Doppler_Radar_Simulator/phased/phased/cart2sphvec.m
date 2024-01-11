function As = cart2sphvec(Ar,az,el)
%cart2sphvec Convert vector from Cartesian to Spherical representation
%   VS = cart2sphvec(VR,AZ,EL) converts the vector VR from its Cartesian
%   representation to its corresponding spherical representation, VS. 
%
%   The spherical representation of a vector is its projection into a local
%   orthonormal basis defined by increasing values of azimuth angle
%   (az_hat), elevation angle (el_hat), and radius (r_hat). The basis
%   varies according its location on a sphere which is determined by the
%   azimuth and elevation angles, AZ and EL (in degrees). AZ must be
%   between -180 and 180 and EL must be between -90 and 90.
%
%   VR is a 3-row matrix whose columns are vectors in the Cartesian
%   representation, in the form of [x;y;z]. VS has the same dimensions as
%   VR. Each column of VS contains the spherical representation, in the
%   form of [az_hat;el_hat;r_hat], of the corresponding column of VR.
%
%   % Example:
%   %   Determine the spherical representation of a vector, [1;0;0],
%   %   defined in the Cartesian coordinate. The spherical coordinate is at
%   %   30 degrees azimuth and 0 degrees elevation.
%
%   vs = cart2sphvec([1;0;0],30,0)

%   Copyright 2012 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

%   translates a vector in rectangular coordinates to its corresponding
%   spherical representation in the coordinates (az_hat, el_hat, r_hat )
%   defined in the direction of (az;el).

% Ar is 3xN rectangular coordinates [x;y;z]
% As is 3xN spherical coordinates [az;el;r] (see cart2sph)
% az is x->y in degrees, scalar
% el is xy->z in degrees, scalar

eml_assert_no_varsize(2:3,Ar,az,el);
validateattributes(Ar,{'double'},{'finite','nonnan','nonempty',...
    '2d','nrows',3},'cart2sphvec','VR');
sigdatatypes.validateAngle(az,'cart2sphvec','AZ',...
    {'scalar','<=',180,'>=',-180});
sigdatatypes.validateAngle(el,'cart2sphvec','EL',...
    {'scalar','<=',90,'>=',-90});
M = phased.internal.azelcoord(az,el);

As = M.'*Ar;

% [EOF]
