function gCoord = local2globalcoord(lclCoord,option,localOrigin,localAxes)
%This function is for internal use only. It may be removed in the future.

%local2globalcoord Local to global coordinates conversion
%   GCOORD = phased.internal.local2globalcoord(LCOORD,OPTION,LORIGIN,LAXES)
%   converts the local coordinate LCOORD to its global counterpart GCOORD.
%   The global coordinates' origin is located at [0; 0; 0] and its three
%   axes are in the directions of [1; 0; 0], [0; 1; 0] and [0; 0; 1],
%   respectively. OPTION specifies the form of LCOORD and GCOORD using one
%   of 'rr' | 'rs' | 'ss' | 'sr', where the default is 'rr'. The 'rr'
%   option converts local rectangular coordinates to global rectangular
%   coordinates. The 'rs' option converts local rectangular coordinates to
%   global spherical coordinates. The 'sr' option converts local spherical
%   coordinates to global rectangular coordinates. The 'ss' option converts
%   local spherical coordinates to global spherical coordinates. LCOORD and
%   GCOORD are 3-row matrices where each column is coordinates in either
%   rectangular or spherical coordinate form.
%
%   If the coordinates are in rectangular form, the three elements
%   represent (x,y,z) in meters. If the coordinates are in spherical form,
%   the three elements represent (az, el, r), where azimuth (az, in
%   degrees) is measured from x axis toward y axis, elevation (el, in
%   degrees) is measured from x-y plane toward z axis and r (in meters) is
%   the radius. If LCOORD is a matrix then the resulting GCOORD is also a
%   matrix with the same dimensions. Each row of GCOORD contains the three
%   coordinates that defines the converted point.
%
%   LORIGIN specifies the origin of the local coordinate system. LORIGIN is
%   a 3x1 column vector containing the rectangular coordinates of the local
%   coordinate system origin with respect to the global coordinate system.
%   The default for LORIGIN is [0; 0; 0];
%
%   LAXES specifies the axes of the local coordinate system. LAXES is a 3x3
%   matrices with each row specifying the direction of local x, y and z
%   axis in the global coordinate system. The default for LAXES is [1 0 0;0
%   1 0;0 0 1]
%
%   % Example:
%   %   A point is located at (2,3,0) in a local coordinate system whose
%   %   origin is at (1,1,0) and whose axes are [0; 1; 0], [1; 0; 0] and 
%   %   [0; 0; 1]. Determine the coordinates of this point in the global 
%   %   system.
%
%   gcoord = phased.internal.local2globalcoord([2; 3; 0],'rr',...
%               [1; 1; 0],[0 1 0;1 0 0;0 0 1])
%
%   See also phased, local2globalcoord.

%   Copyright 2010-2018 The MathWorks, Inc.

%#codegen 

option = convertStringsToChars(option);

if option(1) == 'r'
    lclRect = lclCoord;
else
    lclRect = zeros(size(lclCoord));
    lclCoord(1:2,:) = phased.internal.deg2rad(lclCoord(1:2,:));
    [lclRect(1,:), lclRect(2,:), lclRect(3,:)] = sph2cart(...
        lclCoord(1,:), lclCoord(2,:), lclCoord(3,:));
end

% rotate so that the system is aligned with the global system the local [1;
% 0; 0] should map to [ux; uy; uz] in global note localAxes is [ux vx wx;uy
% vy wy;uz vz wz]; add displacement incurred by the platform position

gRect = localOrigin*ones(1,size(lclRect,2)) + localAxes * lclRect;

if option(2) == 'r' 
    gCoord = gRect;
else
    gCoord = zeros(size(gRect));
    [gCoord(1,:), gCoord(2,:), gCoord(3,:)] = cart2sph(...
        gRect(1,:), gRect(2,:), gRect(3,:));
    gCoord(1:2,:) = phased.internal.rad2deg(gCoord(1:2,:));
end

% [EOF]
