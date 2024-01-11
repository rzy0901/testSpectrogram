function [rcs,az,el] = helperCylinderRCSPattern(c,fc,az,el)
% This function helperCylinderRCSPattern is only in support of
% TargetRCSExample and slexWidebandMonostaticRadarExample. It may change in
% a future release.

% Copyright 2015-2016 The MathWorks, Inc.

% Reference
% [1] Bassem Mahafza, Radar Systems Analysis and Design Using MATLAB, 2nd
% Ed. Chapman & Hall/CRC, 2005

if nargin<3 || isempty(az)
    az = -180:180;
end
if nargin<4 || isempty(el)
    el = -90:90;
end

r = 1;  % cylinder radius
H = 10;  % cylinder height

[elg,~,fcg] = ndgrid(el(:),az(:),fc(:));
lambda = c./fcg;

rcs = lambda.*(r*cosd(elg))./(8*pi*sind(elg).^2);
ind0 = elg==0;
rcs(ind0) = 2*pi*H^2*r./lambda(ind0); % [1] eq 13.50 & 13.51
