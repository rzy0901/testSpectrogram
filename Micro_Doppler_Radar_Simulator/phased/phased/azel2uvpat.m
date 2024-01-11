function [pat_new,ang1_new,ang2_new] = azel2uvpat(pat_old,ang1_old,ang2_old,varargin)
%AZEL2UVPAT Convert radiation pattern from az/el to u/v format
%   PAT_UV = azel2uvpat(PAT_AZEL,AZ,EL) converts the radiation pattern
%   PAT_AZEL represented in azimuth/elevation space to u/v space. AZ and EL
%   are length-P and length-Q vectors, respectively, representing the
%   azimuth and elevation angles (in degrees) at which the original pattern
%   PAT_AZEL is sampled. AZ must be between -90 and 90, and EL must be
%   between -90 and 90. PAT_AZEL is a QxP matrix. 
%   
%   PAT_UV is the converted pattern in u/v space, covering u values from
%   -1 to 1 and v values from -1 to 1. The pattern is uniformly sampled
%   with a step size of 0.01 for both u and v. Note that if the point
%   represented by a pair of U and V is outside the unit circle in u/v
%   space, the corresponding value of PAT_UV is NaN.
%
%   PAT_UV = azel2uvpat(...,U,V) specifies the grid at which PAT_UV is
%   sampled in u/v space. U and V are length-L and length-M vectors,
%   respectively. Both U and V must be between -1 and 1. PAT_UV is an MxL
%   matrix.
%
%   [PAT_UV,U,V] = azel2uvpat(...) also returns the grid U and V in u/v
%   space.
%
%   % Example:
%   %   Convert a radiation pattern in az/el space to u/v space.
%
%   az = -90:90; el = -90:90; pat_azel = ones(181,181);
%   [pat_uv,u,v] = azel2uvpat(pat_azel,az,el);
%   surf(u,v,pat_uv); xlabel('u'); ylabel('v'); zlabel('Pattern');
%
%   See also azel2uv, uv2azel, uv2azelpat, azel2phithetapat.

%   Copyright 2011 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,5,nargin);
eml_assert_no_varsize(1:nargin, pat_old,ang1_old,ang2_old,varargin{:});

if nargin > 4
    ang2_new = varargin{2};
    validateattributes(ang2_new,{'double'},{'vector','real','>=',-1,'<=',1},...
        'azel2uvpat','V');
else
    ang2_new = -1:0.01:1;
end
if nargin > 3
    ang1_new = varargin{1};
    validateattributes(ang1_new,{'double'},{'vector','real','>=',-1,'<=',1},...
        'azel2uvpat','U');
else
    ang1_new = -1:0.01:1;
end

validateattributes(ang1_old,{'double'},{'vector','real','>=',-90,'<=',90},...
    'azel2uvpat','AZ');
validateattributes(ang2_old,{'double'},{'vector','real','>=',-90,'<=',90},...
    'azel2uvpat','EL');
validateattributes(pat_old,{'double'},{'2d','real','size',[numel(ang2_old) numel(ang1_old)]},...
    'azel2uvpat','PAT_AZEL');

fh = @uv2azel;
pat_new = phased.internal.convertpat(fh,pat_old,ang1_old(:),ang2_old(:),ang1_new(:),ang2_new(:),true);




% [EOF]
