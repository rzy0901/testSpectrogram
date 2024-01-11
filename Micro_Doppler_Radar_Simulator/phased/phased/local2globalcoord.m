function gCoord = local2globalcoord(lclCoord,optionArg,localOriginArg,localAxesArg)
%local2globalcoord Local to global coordinates conversion
%   GCOORD = local2globalcoord(LCOORD, OPTION) converts the local
%   coordinate LCOORD to its global counterpart GCOORD.  The global
%   coordinates' origin is located at [0; 0; 0] and its three axes are in
%   the directions of [1; 0; 0], [0; 1; 0] and [0; 0; 1], respectively.
%   OPTION specifies the form of LCOORD and GCOORD using one of 'rr' | 'rs'
%   | 'ss' | 'sr', where the default is 'rr'. The 'rr' option converts
%   local rectangular coordinates to global rectangular coordinates. The
%   'rs' option converts local rectangular coordinates to global spherical
%   coordinates. The 'sr' option converts local spherical coordinates to
%   global rectangular coordinates. The 'ss' option converts local
%   spherical coordinates to global spherical coordinates. LCOORD and
%   GCOORD are 3-row matrices where each column is coordinates in either
%   rectangular or spherical coordinate form.
%
%   If the coordinates are in rectangular form, the three elements
%   represent (x,y,z) in meters. If the coordinates are in spherical form,
%   the three elements represent (az, el, r), where azimuth (az, in
%   degrees) is measured from x axis toward y axis, elevation (el, in
%   degrees) is measured from x-y plane toward z axis and r (in meters) is
%   the radius. If LCOORD is a matrix then the resulting GCOORD is also a
%   matrix with the same dimensions. Each column of GCOORD contains the
%   three coordinates that define the converted point.
%
%   GCOORD = local2globalcoord(...,LORIGIN) specifies the origin, LORIGIN,
%   of the local coordinate system. LORIGIN is a 3x1 column vector
%   containing the rectangular coordinates of the local coordinate system
%   origin with respect to the global coordinate system. The default for
%   LORIGIN is [0; 0; 0].
%
%   GCOORD = local2globalcoord(...,LAXES) specifies the axes, in LAXES, of
%   the local coordinate system. LAXES is a 3x3 matrix with each column
%   specifying the direction of local x, y and z axis in the global
%   coordinate system. The default for LAXES is [1 0 0;0 1 0;0 0 1].
%
%   % Example:
%   %   A point is located at (2,3,0) in a local coordinate system whose
%   %   origin is at (1,1,0) and whose axes are [0; 1; 0], [1; 0; 0] and 
%   %   [0; 0; -1]. Determine the coordinates of this point in the global 
%   %   system.
%
%   gcoord = local2globalcoord([2; 3; 0],'rr',[1; 1; 0],...
%               [0 1 0;1 0 0;0 0 -1])
%
%   See also phased, global2localcoord, rangeangle.

%   Copyright 2008-2011 The MathWorks, Inc.

%   Reference
%   [1] J. D. Foley et al. Computer Graphics: Principles and Practice in C,
%       2nd Ed., Addison-Wesley, 1995

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(1,4,nargin);

if nargin < 4 || isempty(localAxesArg)
    localAxes = eye(3);
else
    localAxes = localAxesArg;
end
if nargin < 3 || isempty(localOriginArg)
    localOrigin = [0; 0; 0];
else
    localOrigin = localOriginArg;
end
if nargin < 2 || isempty(optionArg)
    option = 'rr';
else
    option = optionArg;
end
eml_assert_no_varsize(2:nargin, lclCoord,option,localOrigin,localAxes);
option = validatestring(option,{'rs','rr','sr','ss'},...
    'local2globalcoord','OPTION');
sigdatatypes.validate3DCartCoord(lclCoord,'local2globalcoord','LCOORD');
sigdatatypes.validate3DCartCoord(localOrigin,'local2globalcoord',...
    'LORIGIN',{'size',[3 1]});
sigdatatypes.validate3DCartCoord(localAxes,'local2globalcoord','LAXES',...
    {'size',[3 3]});

for m = 1:3
    localAxes(:,m) = localAxes(:,m)/norm(localAxes(:,m));
end

gCoord = phased.internal.local2globalcoord(lclCoord,option,...
    localOrigin,localAxes);



% [EOF]
