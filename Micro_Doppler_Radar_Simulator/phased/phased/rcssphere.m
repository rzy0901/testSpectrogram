function [rcs,az,el] = rcssphere(r,c,fc,az,el)
%rcssphere   Radar cross section for sphere
%   RCSPAT = rcssphere(R,C,FC) returns the radar cross section (RCS) for a
%   sphere at given frequencies. C is a scalar representing the propagation
%   speed (in m/s) and FC is an L-element vector containing operating
%   frequencies (in Hz) at which the RCS patterns are computed.
%
%   The sphere is located at the origin with a radius of R (in meters).
%
%   RCSPAT = rcssphere(...,AZ,EL) specifies the azimuth angles AZ (in
%   degrees) as a Q-element array and the elevation angles EL (in degrees)
%   as a P-element array. The pairs between azimuth and elevation angles
%   define the directions at which the RCS patterns are computed. The
%   default value for AZ is -180:180 and the default value for EL is
%   -90:90.
%
%   RCSPAT is a PxQxL array containing the RCS pattern (in square meters)
%   of the sphere at given azimuth angles in AZ, elevation angles in EL,
%   and frequencies in FC.
%
%   [RCSPAT,AZ,EL] = rcssphere(...) also returns the azimuth and elevation
%   angles at which the RCS patterns are computed.
%
%   % Example:
%   %   Plot how the RCS of a perfectly conducting sphere vary with
%   %   frequencies. Assume the radius is 10 cm.
%
%   c = 3e8;
%   r = 0.1;
%   kr = 0.3:0.05:30;
%   fc = c*kr/r/(2*pi);
%   rcspat = rcssphere(r,c,fc,0,0);
%   plot(kr,pow2db(rcspat(:)/(pi*r^2)));
%   xlabel('2\pi r/\lambda');
%   ylabel('\sigma/(\pi r^2) (dB)');
%   title('Sphere RCS'); grid on;
%
%   See also phased, rcscylinder, rcstruncone.

% Copyright 2017-2019 The MathWorks, Inc.

% Reference
% [1] Lamont Blake, Calculation of the Radar Cross Section of a Perfectly
% Conducting Sphere, NRL Memorandum Report 2419, 1972
% 
% [2] Roger J. Sullivan, Microwave Radar: Imaging and Advanced Concepts,
% Artech House, 2000

narginchk(3,5);

sigdatatypes.validateDistance(r,'rcssphere','R',{'scalar','positive'});
sigdatatypes.validateSpeed(c,'rcssphere','C',{'scalar','positive'});
sigdatatypes.validateFrequency(fc,'rcssphere','FC',{'vector','positive'});

if nargin<4 || isempty(az)
    az = -180:180;
else
    sigdatatypes.validateAngle(az,'rcssphere','AZ',{'vector','>=',-180,'<=',180});
end
if nargin<5 || isempty(el)
    el = -90:90;
else
    sigdatatypes.validateAngle(el,'rcssphere','EL',{'vector','>=',-90,'<=',90});
end

% A sphere RCS return has three regions, Rayleigh, Mie, and Optical.
% Normally Rayleigh is when the radius is less than a wavelength and
% optical is when the radius is larger than 10 wavelengths. 

% This function simplifies it and groups Mie and Optical together

[~,~,fcg] = ndgrid(el(:),az(:),fc(:));

% since sphere RCS only depends on lambda, each frequency shares the same
% RCS omnidirectional

lambda = c./fc(:).';
k = 2*pi./lambda;
kr = k*r;   % r is scalar

errthreshold = 1e-6;
nrcs_vec = zeros(1,numel(fc));
err = ones(numel(nrcs_vec),1);

n = 0;
idx = (err > errthreshold);
while any(idx)
    n = n+1;
    nrcs_pre = nrcs_vec;
    
    % [1] equation (1)
    nrcs_vec(idx) = nrcs_vec(idx)+(-1)^n*(2*n+1)*(...
        besselj(n+0.5,kr(idx))./besselh(n+0.5,2,kr(idx))+...  %an
        (-(besselj(n+0.5,kr(idx))+kr(idx).*(n*besselj(n-1+0.5,kr(idx))-(n+1)*besselj(n+1+0.5,kr(idx)))./(2*n+1))./...  %bn
        (besselh(n+0.5,2,kr(idx))+kr(idx).*(n*besselh(n-1+0.5,2,kr(idx))-(n+1)*besselh(n+1+0.5,2,kr(idx)))./(2*n+1))));
    err = abs(nrcs_vec(:)-nrcs_pre(:));
    idx = (err > errthreshold);
end

nrcs_vec = abs(nrcs_vec).^2./(kr.^2);
rcs = ones(size(fcg)).*pi*r.^2.*permute(nrcs_vec,[1 3 2]);


