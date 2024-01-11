function [rcs,az,el] = rcsdisc(r,c,fc,az,el)
%rcsdisc   Radar cross section for circular plate
%   RCSPAT = rcsdisc(R,C,FC) returns the radar cross section (RCS) for
%   a circular plate at given frequencies. C is a scalar representing the
%   propagation speed (in m/s) and FC is an L-element vector containing
%   operating frequencies (in Hz) at which the RCS patterns are computed.
%
%   The plate lies on the xy-plane and is centered at the origin with a
%   radius of R (in meters).
%
%   RCSPAT = rcsdisc(...,AZ,EL) specifies the azimuth angles AZ (in
%   degrees) as a Q-element array and the elevation angles EL (in degrees)
%   as a P-element array. The pairs between azimuth and elevation angles
%   define the directions at which the RCS patterns are computed. The
%   default value for AZ is -180:180 and the default value for EL is
%   -90:90.
%
%   RCSPAT is a PxQxL array containing the RCS pattern of the circular
%   plate at given azimuth angles in AZ, elevation angles in EL, and
%   frequencies in FC.
%
%   [RCSPAT,AZ,EL] = rcsdisc(...) also returns the azimuth and
%   elevation angles at which the RCS patterns are computed.
%
%   % Example:
%   %   Plot how the RCS of a circular plate vary along elevation direction
%   %   at 10 GHz. Assume the radius is 40 cm.
%
%   c = 3e8;
%   fc = 10e9;
%   r = 0.4;
%   el = -90:90;
%   rcspat = rcsdisc(r,c,fc,0,el);
%   plot(el,pow2db(rcspat(:)));
%   xlabel('Elevation Angle (GHz)');
%   ylabel('RCS (dBsm)');
%   title('Circular Plate RCS'); grid on;
%
%   See also phased, rcscylinder, rcstruncone.

% Copyright 2017-2018 The MathWorks, Inc.

% Reference
% [1] Bassem Mahafza, Radar Systems Analysis and Design Using MATLAB, 2nd
% Ed. Chapman & Hall/CRC, 2005
%
% [2] Eugene Knott, Radar Cross Section Measurements, SciTech Publishing,
% 2006

%#codegen

narginchk(3,5);

sigdatatypes.validateDistance(r,'rcsdisc','R',{'scalar','positive'});
sigdatatypes.validateSpeed(c,'rcsdisc','C',{'scalar','positive'});
sigdatatypes.validateFrequency(fc,'rcsdisc','FC',{'vector','positive'});

if nargin<4 || isempty(az)
    az = -180:180;
else
    sigdatatypes.validateAngle(az,'rcsdisc','AZ',{'vector','>=',-180,'<=',180});
end
if nargin<5 || isempty(el)
    el = -90:90;
else
    sigdatatypes.validateAngle(el,'rcsdisc','EL',{'vector','>=',-90,'<=',90});
end

[elg,azg,fcg] = ndgrid(el(:),az(:),fc(:));
lambda = c./fcg;
azg = deg2rad(azg); %#ok<NASGU>
elg = deg2rad(elg);

kr = 2*pi*r./lambda;

% [1] describes the equation as spherical Bessel function, which is
% beseelj(1.5,...), but [2] seems to suggest it's just regular Bessel
% function. The code in [1] also uses regular Bessel function, so we use
% regular Bessel function here.
rcs = pi*kr.^2.*r^2.*(2*besselj(1,2*kr.*cos(elg))./(2*kr.*cos(elg))).^2.*sin(elg).^2;  % [1] eq 13.39
ind0 = (elg==-pi/2|elg==pi/2);
rcs(ind0) = 4*pi^3*r^4./lambda(ind0).^2;  % [1] eq 13.37
