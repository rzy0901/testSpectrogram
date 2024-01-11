function [pat_new,ang1_new,ang2_new] = azel2phithetapat(pat_old,ang1_old,ang2_old,varargin)
%AZEL2PHITHETAPAT Convert radiation pattern from az/el to phi/theta format
%   PAT_PHITHETA = azel2phithetapat(PAT_AZEL,AZ,EL) converts the radiation
%   pattern PAT_AZEL represented in azimuth/elevation space to phi/theta
%   space. AZ and EL are length-P and length-Q vectors, respectively,
%   representing the az and el values at which the original pattern
%   PAT_AZEL is sampled. AZ must be between -180 and 180, and EL must be
%   between -90 and 90. PAT_AZEL is a QxP matrix.
%
%   PAT_PHITHETA is the converted pattern in phi/theta space, covering phi
%   values from 0 to 360 and theta values from 0 to 180 (in degrees). The
%   pattern is uniformly sampled with a step size of 1 in both phi and
%   theta.
%
%   PAT_PHITHETA = azel2phithetapat(...,PHI,THETA) specifies the grid at
%   which PAT_PHITHETA is sampled in phi/theta space (in degrees). PHI
%   and THETA are length-L and length-M vectors, respectively. PHI must be
%   between 0 and 360 and THETA must be between 0 and 180. PAT_PHITHETA
%   is an MxL matrix.
%
%   [PAT_PHITHETA,PHI,THETA] = azel2phithetapat(...) also returns the grid
%   PHI and THETA in phi/theta space.
%
%   % Example:
%   %   Convert a radiation pattern in az/el space to phi/theta space.
%
%   az = -180:180; el = -90:90; pat_azel = ones(181,361);
%   [pat_phitheta,phi,theta] = azel2phithetapat(pat_azel,az,el);
%   surf(phi,theta,pat_phitheta); 
%   xlabel('phi (degrees)'); ylabel('theta (degrees)'); zlabel('Pattern');
%
%   See also azel2phitheta, phitheta2azel, phitheta2azelpat,
%   azel2uvpat.

%   Copyright 2011 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,5,nargin);
eml_assert_no_varsize(1:nargin, pat_old,ang1_old,ang2_old,varargin{:});

if nargin > 4
    ang2_new = varargin{2};
    validateattributes(ang2_new,{'double'},{'vector','real','>=',0,'<=',180},...
        'azel2phithetapat','THETA');
else
    ang2_new = 0:180;
end
if nargin > 3
    ang1_new = varargin{1};
    validateattributes(ang1_new,{'double'},{'vector','real','>=',0,'<=',360},...
        'azel2phithetapat','PHI');
else
    ang1_new = 0:360;
end

validateattributes(ang1_old,{'double'},{'vector','real','>=',-180,'<=',180},...
    'azel2phithetapat','AZ');
validateattributes(ang2_old,{'double'},{'vector','real','>=',-90,'<=',90},...
    'azel2phithetapat','EL');
validateattributes(pat_old,{'double'},{'2d','real','size',[numel(ang2_old) numel(ang1_old)]},...
    'azel2phithetapat','PAT_AZEL');

fh = @phitheta2azel;
pat_new = phased.internal.convertpat(fh,pat_old,ang1_old(:),ang2_old(:),ang1_new(:),ang2_new(:),false);


% [EOF]
