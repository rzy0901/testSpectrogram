function lclCoord = global2localcoord(gCoord,option,localOrigin,localAxes)
%This function is for internal use only. It may be removed in the future.

%global2localcoord Global to local coordinates conversion
%   LCOORD = phased.internal.global2localcoord(GCOORD,OPTION,LORIGIN,LAXES)
%   converts the global coordinate GCOORD to its local counterpart LCOORD.
%   The global coordinates' origin is located at [0; 0; 0] and its three
%   axes are in the directions of [1; 0; 0], [0; 1; 0] and [0; 0; 1],
%   respectively. OPTION specifies the form of GCOORD and LCOORD using one
%   of 'rr' | 'rs' | 'ss' | 'sr', where the default is 'rr'. The 'rr'
%   option converts global rectangular coordinates to local rectangular
%   coordinates. The 'rs' option converts global rectangular coordinates to
%   local spherical coordinates. The 'sr' option converts global spherical
%   coordinates to local rectangular coordinates. The 'ss' option converts
%   global spherical coordinates to local spherical coordinates. GCOORD and
%   LCOORD are 3 row matrices where each column is coordinates in either
%   rectangular or spherical coordinate form.
%
%   If the coordinates are in rectangular form, the three elements
%   represent (x,y,z) in meters. If the coordinates are in spherical form,
%   the three elements represent (az, el, r), where azimuth (az, in
%   degrees) is measured from x axis toward y axis, elevation (el, in
%   degrees) is measured from x-y plane toward z axis and r (in meters) is
%   the radius. If GCOORD is a matrix, the resulting LCOORD is also a
%   matrix with the same dimensions. Each column of LCOORD contains the
%   three coordinates that defines the converted point.
%
%   LORIGIN specifies the origin of the local coordinate system. LORIGIN is
%   a 3x1 column vector containing the rectangular coordinates of the local
%   coordinate system origin with respect to the global coordinate system.
%   The default for LORIGIN is [0; 0; 0].
%
%   LAXES specifies the axes of the local coordinate system. LAXES is a 3x3
%   matrices with each row specifying the direction of local x, y and z
%   axis in the global coordinate system. The default for LAXES is [1 0 0;0
%   1 0;0 0 1].
%
%   % Example:
%   %   A point is located at (2,3,0) in the global coordinate system.
%   %   Determine the coordinates of this point in the local coordinate 
%   %   system whose origin is at (1,1,0) and whose axes are [0; 1; 0], 
%   %   [1; 0; 0] and [0; 0; 1].
%
%   lclcoord = phased.internal.global2localcoord([2; 3; 0],'rr',...
%               [1; 1; 0],[0 1 0;1 0 0;0 0 1])
%
%   See also phased, global2localcoord.

%   Copyright 2010-2018 The MathWorks, Inc.

% internal function, no error checking is performed

%#codegen 

option = convertStringsToChars(option);

if option(1) == 'r' 
    gRect = gCoord;
else
    gRect = zeros(size(gCoord));
    gCoord(1:2,:) = phased.internal.deg2rad(gCoord(1:2,:));
    [gRect(1,:), gRect(2,:), gRect(3,:)] = sph2cart(...
        gCoord(1,:), gCoord(2,:), gCoord(3,:));
end

% subtract displacement incurred by the platform position; rotate so that
% the system is local the global [ux; uy; uz] should map to [1; 0; 0] in local
% note localAxes is [ux vx wx;uy vy wy;uz vz wz];

lclRect = localAxes.' * (gRect - localOrigin*ones(1,size(gRect,2)));

if option(2) == 'r' 
    lclCoord = lclRect;
else
    lclCoord = zeros(size(lclRect));
    [lclCoord(1,:), lclCoord(2,:), lclCoord(3,:)] = cart2sph(...
        lclRect(1,:), lclRect(2,:), lclRect(3,:));
    lclCoord(1:2,:) = phased.internal.rad2deg(lclCoord(1:2,:));
end


% [EOF]
