function [rcs,az,el] = rcstruncone(r1,r2,H,c,fc,az,el)
%rcstruncone   Radar cross section for frustum
%   RCSPAT = rcstruncone(R1,R2,H,C,FC) returns the radar cross section
%   (RCS) for a frustum (truncated cone) at given frequencies. C is a
%   scalar representing the propagation speed (in m/s) and FC is an
%   L-element vector containing operating frequencies (in Hz) at which the
%   RCS patterns are computed.
%
%   The frustum is assumed to be along the z axis. R1 is the radius of the
%   small end of the frustum and r2 is the radius of the big end of the
%   frustum. H is the height of the frustum. The origin is located at the
%   tip of the cone if it were not cut. Note that the frustum is not
%   covered in the two ends so there is no return for elevation angles
%   close to -90 and 90 degrees.
%
%   RCSPAT = rcstruncone(...,AZ,EL) specifies the azimuth angles AZ (in
%   degrees) as a Q-element array and the elevation angles EL (in degrees)
%   as a P-element array. The pairs between azimuth and elevation angles
%   define the directions at which the RCS patterns are computed. The
%   default value for AZ is -180:180 and the default value for EL is
%   -90:90.
%
%   RCSPAT is a PxQxL array containing the RCS pattern of the frustum at
%   given azimuth angles in AZ, elevation angles in EL, and frequencies in
%   FC.
%
%   [RCSPAT,AZ,EL] = rcstruncone(...) also returns the azimuth and
%   elevation angles at which the RCS patterns are computed.
%
%   % Example:
%   %   Compute the RCS pattern for a frustum whose small end has a radius
%   %   of 2 cm, large end has a radius of 6 cm, and a height of 20 cm. 
%   %   Assume the operating frequency is 9 GHz. Plot the pattern along 
%   %   elevation.
%
%   c = 3e8;
%   fc = 9e9;
%   r1 = 0.02;
%   r2 = 0.06;
%   h = 0.2;
%   el = -89.5:0.5:89.5;
%   rcspat = rcstruncone(r1,r2,h,c,fc,0,el);
%   plot(el,pow2db(rcspat));
%   xlabel('Elevation Angle (deg)');
%   ylabel('RCS (dBsm)');
%   title('Frustum RCS'); grid on;
%
%   See also phased, rcssphere, rcscylinder.

% Copyright 2017-2018 The MathWorks, Inc.

% Reference
% [1] Bassem Mahafza, Radar Systems Analysis and Design Using MATLAB, 2nd
% Ed. Chapman & Hall/CRC, 2005

%#codegen

narginchk(5,7);

sigdatatypes.validateDistance(r1,'rcstruncone','R1',{'scalar'}); % can be 0 for cone
sigdatatypes.validateDistance(r2,'rcstruncone','R2',{'scalar','positive','>',r1});
sigdatatypes.validateDistance(H,'rcstruncone','H',{'scalar','positive'});
sigdatatypes.validateSpeed(c,'rcstruncone','C',{'scalar','positive'});
sigdatatypes.validateFrequency(fc,'rcstruncone','FC',{'vector','positive'});

if nargin<6 || isempty(az)
    az = -180:180;
else
    sigdatatypes.validateAngle(az,'rcstruncone','AZ',{'vector','>=',-180,'<=',180});
end
if nargin<7 || isempty(el)
    el = -90:90;
else
    sigdatatypes.validateAngle(el,'rcstruncone','EL',{'vector','>=',-90,'<=',90});
end

[elg,azg,fcg] = ndgrid(el(:),az(:),fc(:));
lambda = c./fcg;
azg = deg2rad(azg); %#ok<NASGU>
elg = deg2rad(elg);

tanalpha = (r2-r1)/H;
z1 = r1/tanalpha;
z2 = z1+H;

% non-normal incidence
alpha = atan(tanalpha);
rcs = lambda.*tanalpha./(8*pi*sin(elg+pi/2)).*tan(elg+pi/2+alpha).^2; % [1] eq 13.47
smallendind = elg<-alpha;
rcs(smallendind) = rcs(smallendind)*z1;
bigendind = elg>-alpha;              
rcs(bigendind) = rcs(bigendind)*z2;

% normal incidence
ind0 = abs(elg-(-alpha))<eps;
rcs(ind0) = 8*pi*(z2^(3/2)-z1^(3/2)).^2./(9*lambda(ind0).*sin(elg(ind0)+pi/2)).*tanalpha.*...
    (sin(elg(ind0)+pi/2)+cos(elg(ind0)+pi/2).*tanalpha).^2; % [1] eq 13.44 

% two ends
tol = 0.001;
inde1 = elg<deg2rad(-90+tol);
rcs(inde1) = 0;
inde2 = elg>deg2rad(90-tol);
rcs(inde2) = 0;
