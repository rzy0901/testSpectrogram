function [angUpperL,angLowerL,angUpperR,angLowerR] = getBicyclistLegAngles(spacing,mounts,upperLegVec,lowerLegVec)
%This function is for internal use only. It may be removed in the future.

%   Copyright 2019 The MathWorks, Inc.

%#codegen


%% Left leg
% Intersect Circles
distHipToFootLeft = norm(mounts.Hip(:,1)-mounts.Feet(:,1));
R = distHipToFootLeft;
x1 = mounts.Hip(1,1);
z1 = mounts.Hip(3,1);
r1 = upperLegVec(end);
x2 = mounts.Feet(1,1);
z2 = mounts.Feet(3,1);
r2 = lowerLegVec(end)+spacing;
[xLknee,zLknee] = intersectCircles(R,r1,x1,z1,r2,x2,z2);
angUpperL = atan((mounts.Hip(3,1)-zLknee)/(mounts.Hip(1,1)-xLknee));
angLowerL = atan((zLknee-mounts.Feet(3,1))/(xLknee-mounts.Feet(1,1)));


%% Right leg
% Intersect Circles
distHipToFootRight = norm(mounts.Hip(:,2)-mounts.Feet(:,2));
R = distHipToFootRight;
x1 = mounts.Hip(1,2);
z1 = mounts.Hip(3,2);
r1 = upperLegVec(end);
x2 = mounts.Feet(1,2);
z2 = mounts.Feet(3,2);
r2 = lowerLegVec(end)+spacing;
[xRknee,zRknee] = intersectCircles(R,r1,x1,z1,r2,x2,z2);
angUpperR = atan((mounts.Hip(3,2)-zRknee)/(mounts.Hip(1,2)-xRknee));
angLowerR = atan((zRknee-mounts.Feet(3,2))/(xRknee-mounts.Feet(1,2)));
end


%% Supporting functions
function [x,y] = intersectCircles(R,r1,x1,y1,r2,x2,y2)
% Find intersection of two circles to find location of knee
a = 1/2;
b = (r1^2-r2^2)/(2*R^2);
c = 1/2*sqrt(2*(r1^2+r2^2)/(R^2)-(r1^2-r2^2)^2/(R^4)-1);
x = a*(x1+x2)+b*(x2-x1)-c*(y2-y1);
y = a*(y1+y2)+b*(y2-y1)-c*(x1-x2);
end

