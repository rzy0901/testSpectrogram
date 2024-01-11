function [pat_new,ang1_new,ang2_new] = phitheta2azelpat(pat_old,ang1_old,ang2_old,varargin)
%PHITHETA2AZELPAT Convert radiation pattern from phi/theta to az/el format
%   PAT_AZEL = phitheta2azelpat(PAT_PHITHETA,PHI,THETA) converts the
%   radiation pattern PAT_PHITHETA represented in phi/theta space to
%   azimuth/elevation space. PHI and THETA are length-P and length-Q
%   vectors, respectively, representing the phi and theta angles (in
%   degrees) at which the original pattern PAT_PHITHETA is sampled. PHI
%   must be between 0 and 360, and THETA must be between 0 and 180.
%   PAT_PHITHETA is a QxP matrix.
%   
%   PAT_AZEL is the converted pattern in az/el space covering az from -180
%   to 180 and el from -90 to 90 (in degrees). The pattern is uniformly
%   sampled with a step size of 1 for both az and el.
%
%   PAT_AZEL = phitheta2azelpat(...,AZ,EL) specifies the grid at which
%   PAT_AZEL is sampled in az/el space (in degrees). AZ and EL are
%   length-L and length-M vectors, respectively. AZ must be between -180
%   and 180 and EL must be between -90 and 90. PAT_AZEL is an MxL matrix.
%
%   [PAT_AZEL,AZ,EL] = phitheta2azelpat(...) also returns the grid AZ and
%   EL in az/el space.
%
%   % Example:
%   %   Convert a radiation pattern in phi/theta space to az/el space.
%
%   phi = 0:360; theta = 0:90; pat_phitheta = ones(91,361);
%   [pat_azel,az,el] = phitheta2azelpat(pat_phitheta,phi,theta);
%   surf(az,el,pat_azel); 
%   xlabel('az (degrees)'); ylabel('el (degrees)'); zlabel('Pattern');
%
%   See also phitheta2azel, azel2phitheta, azel2phithetapat,
%   phitheta2uvpat.

%   Copyright 2011 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,5,nargin);
eml_assert_no_varsize(1:nargin,pat_old,ang1_old,ang2_old,varargin{:});
if nargin > 4
    ang2_new = varargin{2};
    validateattributes(ang2_new,{'double'},{'vector','real','>=',-90,'<=',90},...
        'phitheta2azelpat','EL');
else
    ang2_new = -90:90;
end
if nargin > 3
    ang1_new = varargin{1};
    validateattributes(ang1_new,{'double'},{'vector','real','>=',-180,'<=',180},...
        'phitheta2azelpat','AZ');
else
    ang1_new = -180:180;
end

validateattributes(ang1_old,{'double'},{'vector','real','>=',0,'<=',360},...
    'phitheta2azelpat','PHI');
validateattributes(ang2_old,{'double'},{'vector','real','>=',0,'<=',180},...
    'phitheta2azelpat','THETA');
validateattributes(pat_old,{'double'},{'2d','real','size',[numel(ang2_old) numel(ang1_old)]},...
    'phitheta2azelpat','PAT_PHITHETA');

fh = @azel2phitheta;
pat_new = phased.internal.convertpat(fh,pat_old,ang1_old(:),ang2_old(:),ang1_new(:),ang2_new(:),false);


% [EOF]
