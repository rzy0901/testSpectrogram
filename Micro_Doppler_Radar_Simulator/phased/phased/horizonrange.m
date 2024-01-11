function Rh = horizonrange(platformHeight, Re)
%horizonrange Horizon range
%   RH = horizonrange(H) returns the horizon range (in meters), RH, of a
%   radar system located at a height of H (in meters) above the surface.
%   The calculation uses 4/3 of earth radius as the effective earth radius.
%
%   RH = horizonrange(H,RE) specifies the effective earth radius (in
%   meters) as a positive scalar, RE, where the default is 4/3 of the earth
%   radius.
%
%   % Example:
%   %   Determine the horizon range of an antenna that is 30 meters above
%   %   the ground.
%
%   Rh = horizonrange(30)
%
%   See also effearthradius, grazingang, depressionang.

%   Copyright 2008-2011 The MathWorks, Inc.
%     

%   Reference
%   [1] James Ward, Space-Time Adaptive Processing for Airborne Radar, 1994
%   [2] Long, Radar Reflectivity of Land and Sea, 2001

%#codegen 
%#ok<*EMCA>


phased.internal.narginchk(1,2,nargin);

if nargin < 2
    Re = effearthradius;
end
    
sigdatatypes.validateDistance(platformHeight,'horizonrange','H',{'vector'});
sigdatatypes.validateDistance(Re,'horizonrange','RE',{'scalar','positive'});

Rh = sqrt(2.*platformHeight.*Re+platformHeight.^2);



% [EOF]
