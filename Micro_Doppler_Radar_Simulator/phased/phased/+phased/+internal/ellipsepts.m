function epts = ellipsepts(focal1,focal2,r,ang)
%This function is for internal use only. It may be removed in the future.

% This function is only in support of ArrayProcessingForMIMOExample. It may
% be removed in a future release.

%   Copyright 2016 The MathWorks, Inc.

% r is the total path length, note r = 2a
% focal 1, [x;y]
% focal 2, [x;y]
% ang, angles where points are, for full ellipse, use -180;180
% epts, points on the ellipse [x;y]

%#codegen

c = norm(focal1-focal2)/2;
elipang = atand((focal2(2)-focal1(2))/(focal2(1)-focal1(1)));
elipcenter = (focal1+focal2)/2;
a = r/2;
b = sqrt(a^2-c^2);
posx = a*cosd(ang);
posy = b*sind(ang);
epts = [posx;posy];
epts = [cosd(elipang) -sind(elipang);sind(elipang) cosd(elipang)]*epts;
epts = epts+elipcenter;

