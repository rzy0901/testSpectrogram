function Re = effearthradius(rgradient)
%effearthradius Effective earth radius
%   Re = effearthradius returns the effective radius of spherical earth (in
%   meters). The effective radius is calculated using a refractivity
%   gradient of -39e-9, which results in approximately 4/3 of the real
%   earth radius.
%
%   Re = effearthradius(RGradient) specifies the refractivity gradient
%   RGradient (in units of 1/meter).
%
%   % Example:
%   %   Calculate the effective earth radius.
%
%   re = effearthradius
%
%   See also horizonrange, depressionang.

%   Copyright 2010-2011 The MathWorks, Inc.

%   Reference
%   [1] James Ward, Space-Time Adaptive Processing for Airborne Radar, 
%       Lincoln Lab Technical Report, 1994
%   [2] Long, Radar Reflectivity of Land and Sea, 3rd Ed. Artech House, 
%       2001
%   [3] Mahafza, Radar Signal Analysis and Processing Using MATLAB, CRC 
%       Press, 2009
%   [4] Skolnik, Introduction to Radar Systems, 3rd Ed. McGraw Hill, 2002

%   If no refraction, the curvature is 1/R, which means that the horizontal
%   beam is always parallel to the earth.  However, with the change of
%   refractivity for each layer, it's no longer true.  The effective earth
%   radius is the radius that the modified curvature is 1/Re

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(0,1,nargin);


if ~nargin
    rgradient = -39e-9;
else
    eml_assert_no_varsize(1,rgradient);
end

validateattributes(rgradient,{'numeric'},{'scalar','real','finite','<=',0},...
    'effearthradius','RGradient');

R =physconst('earthradius');
k = 1/(1+R*rgradient);
Re = k*R;


% [EOF]
