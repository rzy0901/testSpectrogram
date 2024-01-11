function rspeed = radialspeed(tgtpos, tgtvel, refpos, refvel)
%radialspeed Relative radial speed
%   RSPEED = radialspeed(POS,V) returns the radial speed RSPEED (in meters
%   per second) of the given platforms, with position POS and velocity V,
%   relative to a static reference platform at the origin. POS must be a
%   3xN matrix whose columns specify the positions in [x; y; z] format (in
%   meters). V must also be a 3xN matrix whose columns specify the
%   velocities in [x; y; z] format (in meters per second). RSPEED is an Nx1
%   vector with each element representing the radial speed (in meters per
%   second) corresponding to the position and the velocity specified in POS
%   and V.
%
%   RSPEED = radialspeed(POS,V,REFPOS) specifies the reference platform's
%   position REFPOS. REFPOS is a 3x1 vector specifying the reference
%   platform's position in [x; y; z] format (in meters).
%
%   RSPEED = radialspeed(POS,V,REFPOS,REFV) specifies the reference
%   platform's velocity REFV. REFV is a 3x1 vector specifying the reference
%   platform velocity in [x; y; z] format (in meters per second).
%
%   When RSPEED is positive, the platform is approaching the reference
%   platform. When RSPEED is negative, the platform is moving away from the
%   reference platform.
%
%   % Example:
%   %   Calculate the radial speed of a target relative to a stationary
%   %   platform. Assume that the target is moving at a speed of 
%   %   [10; 10; 0] meters per second and is located at [20; 20; 0] meters.
%   %   The reference platform is located at [1; 1; 0].
%
%   rspeed = radialspeed([20; 20; 0],[10; 10; 0],[1; 1; 0])
%
%   See also phased, speed2dop.

%   Copyright 2009-2012 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,4,nargin);

if nargin < 4
    refvel = [0; 0; 0];
end

if nargin < 3
    refpos = [0; 0; 0];
end

eml_assert_no_varsize(1:nargin, tgtpos, tgtvel, refpos, refvel);
sigdatatypes.validate3DCartCoord(tgtpos,'radialspeed','POS');
sigdatatypes.validate3DCartCoord(tgtvel,'radialspeed','V');
sigdatatypes.validate3DCartCoord(refpos,'radialspeed','REFPOS',{'column'});
sigdatatypes.validate3DCartCoord(refvel,'radialspeed','REFV',{'column'});

tgtnum = size(tgtpos,2);
cond = size(tgtvel,2) == tgtnum;
if ~cond
    coder.internal.assert(cond,'phased:radialspeed:PosVelMismatch');
end

tgtdirec = bsxfun(@minus, tgtpos, refpos);
veldirec = bsxfun(@minus, tgtvel, refvel);

%Take the 2-norm of veldirec and tgtdirec
rn = sqrt(sum(tgtdirec.^2));
sn = sqrt(sum(veldirec.^2));

% negative sign to ensure that incoming relative speed is positive
rspeed = -(sum(veldirec.*tgtdirec)./rn);


%now take care of rspeed corner cases:
rspeed(sn==0) = 0;
% when co-located, the radial speed is relative speed, but departing
rnEq0 = (rn==0);
rspeed(rnEq0) = -sn(rnEq0); 

%make a column vector
rspeed = reshape(rspeed, [],1); 


% [EOF]
