function [rcs,az,el] = rcscylinder(r1,r2,H,c,fc,az,el)
%rcscylinder   Radar cross section for cylinder
%   RCSPAT = rcscylinder(R1,R2,H,C,FC) returns the radar cross section
%   (RCS) for an elliptical cylinder at given frequencies. C is a scalar
%   representing the propagation speed (in m/s) and FC is an L-element
%   vector containing operating frequencies (in Hz) at which the RCS
%   patterns are computed.
%
%   The elliptical cylinder is assumed to be sitting on the x-y plane. The
%   major axis, R1 (in meters), lies along the x axis and the minor axis,
%   R2 (in meters), lies along the y axis. H (in meters) is the height of
%   the cylinder along the z axis. For circular cylinders, set R1 and R2 to
%   be the same.
%
%   RCSPAT = rcscylinder(...,AZ,EL) specifies the azimuth angles AZ (in
%   degrees) as a Q-element array and the elevation angles EL (in degrees)
%   as a P-element array. The pairs between azimuth and elevation angles
%   define the directions at which the RCS patterns are computed. The
%   default value for AZ is -180:180 and the default value for EL is
%   -90:90.
%
%   RCSPAT is a PxQxL array containing the RCS pattern of the cylinder at
%   given azimuth angles in AZ, elevation angles in EL, and frequencies in
%   FC.
%
%   [RCSPAT,AZ,EL] = rcscylinder(...) also returns the azimuth and
%   elevation angles at which the RCS patterns are computed.
%
%   % Example:
%   %   Compute the RCS pattern for a circular cylinder whose radius is 
%   %   12.5 cm and the height is 1 m. Assume the operating frequency is 
%   %   3.5 GHz. Plot the pattern along elevation.
%
%   c = 3e8;
%   fc = 3.5e9;
%   r = 0.125;
%   h = 1;
%   el = -90:90;
%   rcspat = rcscylinder(r,r,h,c,fc,0,el);
%   plot(el,pow2db(rcspat));
%   ylim([-50 10])
%   xlabel('Elevation Angle (deg)');
%   ylabel('RCS (dBsm)');
%   title('Cylinder RCS'); grid on;
%
%   See also phased, rcssphere, rcstruncone.

% Copyright 2017-2018 The MathWorks, Inc.

% Reference
% [1] Bassem Mahafza, Radar Systems Analysis and Design Using MATLAB, 2nd
% Ed. Chapman & Hall/CRC, 2005

%#codegen

narginchk(5,7);

sigdatatypes.validateDistance(r1,'rcscylinder','R1',{'scalar','positive'});
sigdatatypes.validateDistance(r2,'rcscylinder','R2',{'scalar','positive'});
sigdatatypes.validateDistance(H,'rcscylinder','H',{'scalar','positive'});
sigdatatypes.validateSpeed(c,'rcscylinder','C',{'scalar','positive'});
sigdatatypes.validateFrequency(fc,'rcscylinder','FC',{'vector','positive'});

if nargin<6 || isempty(az)
    az = -180:180;
else
    sigdatatypes.validateAngle(az,'rcscylinder','AZ',{'vector','>=',-180,'<=',180});
end
if nargin<7 || isempty(el)
    el = -90:90;
else
    sigdatatypes.validateAngle(el,'rcscylinder','EL',{'vector','>=',-90,'<=',90});
end

[elg,azg,fcg] = ndgrid(el(:),az(:),fc(:));
lambda = c./fcg;
azg = deg2rad(azg);
elg = deg2rad(elg);


rcs = (lambda.*r2.^2.*r1.^2.*cos(elg))./(8*pi*sin(elg).^2.*(r1.^2.*cos(azg).^2+r2.^2.*sin(azg).^2).^1.5);
ind0 = elg==0;
rcs(ind0) = 2*pi*H^2*r2.^2*r1.^2./(lambda(ind0).*(r1.^2.*cos(azg(ind0)).^2+r2.^2.*sin(azg(ind0)).^2).^1.5); % [1] eq 13.50 & 13.51
