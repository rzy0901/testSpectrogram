function [pat_new,ang1_new,ang2_new] = uv2azelpat(pat_old,ang1_old,ang2_old,varargin)
%UV2AZELPAT Convert radiation pattern from u/v to az/el format
%   PAT_AZEL = uv2azelpat(PAT_UV,U,V) converts the radiation pattern PAT_UV
%   represented in u/v space to azimuth/elevation space. U and V are
%   length-P and length-Q vectors, respectively, representing the u and v
%   values at which the original pattern PAT_UV is sampled. Both U and V
%   must be between -1 and 1. PAT_UV is a QxP matrix.
%
%   PAT_AZEL is the converted pattern in az/el space covering az from -90
%   to 90 and el from -90 to 90 (in degrees). The pattern is uniformly
%   sampled with a step size of 1 for both az and el.
%
%   PAT_AZEL = uv2azelpat(...,AZ,EL) specifies the grid at which PAT_AZEL
%   is sampled in az/el space (in degrees). AZ and EL are length-L and
%   length-M vectors, respectively. AZ must be between -90 and 90 and EL
%   must be between -90 and 90. PAT_AZEL is an MxL matrix.
%
%   [PAT_AZEL,AZ,EL] = uv2azelpat(...) also returns the grid AZ and
%   EL in az/el space.
%
%   % Example:
%   %   Convert a radiation pattern in u/v space to az/el space.
%
%   u = -1:0.01:1; v = -1:0.01:1; pat_uv = ones(201,201);
%   [pat_azel,az,el] = uv2azelpat(pat_uv,u,v);
%   surf(az,el,pat_azel); 
%   xlabel('az (degrees)'); ylabel('el (degrees)'); zlabel('Pattern');
%
%   See also uv2azel, azel2uv, azel2uvpat, uv2phithetapat.

%   Copyright 2011 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,5,nargin);
eml_assert_no_varsize(1:nargin,pat_old,ang1_old,ang2_old,varargin{:});
if nargin > 4
    ang2_new = varargin{2};
    validateattributes(ang2_new,{'double'},{'vector','real','>=',-90,'<=',90},...
        'uv2azelpat','EL');
else
    ang2_new = -90:90;
end
if nargin > 3
    ang1_new = varargin{1};
    validateattributes(ang1_new,{'double'},{'vector','real','>=',-90,'<=',90},...
        'uv2azelpat','AZ');
else
    ang1_new = -90:90;
end

validateattributes(ang1_old,{'double'},{'vector','real','>=',-1,'<=',1},...
    'uv2azelpat','U');
validateattributes(ang2_old,{'double'},{'vector','real','>=',-1,'<=',1},...
    'uv2azelpat','V');
validateattributes(pat_old,{'double'},{'2d','real','size',[numel(ang2_old) numel(ang1_old)]},...
    'uv2azelpat','PAT_UV');

fh = @azel2uv;
pat_new = phased.internal.convertpat(fh,pat_old,ang1_old(:),ang2_old(:),ang1_new(:),ang2_new(:),false);


% [EOF]
