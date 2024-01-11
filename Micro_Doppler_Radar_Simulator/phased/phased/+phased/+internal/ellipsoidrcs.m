function [rcs,az,el] = ellipsoidrcs(ra,rb,rc,c,fc,az,el)
%This function is for internal use only. It may be removed in the future.

% Copyright 2017-2018 The MathWorks, Inc.

% Reference
% [1] Bassem Mahafza, Radar Systems Analysis and Design Using MATLAB, 2nd
% Ed. Chapman & Hall/CRC, 2005
%
% [2] David Jenn, Radar and Laser Cross Section Engineering, 2nd Ed.
% American Institute of Aeronautics and Astronautics, 2005

% this result is only valid in optical region where ka >> 1 where k is wave
% number and a is the dimension

%#codegen

narginchk(5,7);

if nargin<6 || isempty(az)
    az = -180:180;
end
if nargin<7 || isempty(el)
    el = -90:90;
end

if isscalar(fc)
    [elg,azg] = ndgrid(el(:),az(:));
    fcg = fc;
else
    [elg,azg,fcg] = ndgrid(el(:),az(:),fc(:));
end
lambda = c./fcg;

k = 10;
ignoreWavelength = true;
% check far field condition
cond = any(max([ra rb rc]) < opticalregion(lambda,k)) && ~ignoreWavelength;
if cond
    fmin = k*c/(2*pi*max([ra rb rc]));
    coder.internal.errorIf(cond,'phased:phased:expectedGreaterThan','FREQ',sprintf('%5.2f',fmin));
end

azg = deg2rad(azg);
elg = deg2rad(elg);

rcs =  pi*ra.^2*rb.^2*rc.^2./...    
    (ra.^2.*cos(elg).^2.*cos(azg).^2+rb.^2.*cos(elg).^2.*sin(azg).^2+rc.^2.*sin(elg).^2).^2;  % [1] eq 13.34

function D = opticalregion(lambda,fac)
%opticalregion    Compute the optical region

% Electrically small is defined as 2*pi*D/lambda < fac
% fac is a factor, default is 10 and shouldn't be less than 5

narginchk(1,2);
if nargin < 2
    fac = 10;
end

D = fac*lambda/(2*pi);
