function [pat_new,ang1_new,ang2_new] = phitheta2uvpat(pat_old,ang1_old,ang2_old,varargin)
%PHITHETA2UVPAT Convert radiation pattern from phi/theta to u/v format
%   PAT_UV = phitheta2uvpat(PAT_PHITHETA,PHI,THETA) converts the radiation
%   pattern PAT_PHITHETA represented in phi/theta space to u/v space. PHI
%   and THETA are length-P and length-Q vectors, respectively, representing
%   the phi and theta angles (in degrees) at which the original pattern
%   PAT_PHITHETA is sampled. PHI must be between 0 and 360, and THETA must
%   be between 0 and 90. PAT_PHITHETA is a QxP matrix.
%   
%   PAT_UV is the converted pattern in u/v space, covering u values from
%   -1 to 1 and v values from -1 to 1. The pattern is uniformly sampled
%   with a step size of 0.01 for both u and v. Note that if the point
%   represented by a pair of U and V is outside the unit circle in u/v
%   space, the corresponding value of PAT_UV is NaN.
%
%   PAT_UV = phitheta2uvpat(...,U,V) specifies the grid at which PAT_UV is
%   sampled in u/v space. U and V are length-L and length-M vectors,
%   respectively. Both U and V must be between -1 and 1. PAT_UV is an MxL
%   matrix.
%
%   [PAT_UV,U,V] = phitheta2uvpat(...) also returns the grid U and V in
%   u/v space.
%
%   % Example:
%   %   Convert a radiation pattern in phi/theta space to u/v space.
%
%   phi = 0:360; theta = 0:90; pat_phitheta = ones(91,361);
%   [pat_uv,u,v] = phitheta2uvpat(pat_phitheta,phi,theta);
%   surf(u,v,pat_uv); xlabel('u'); ylabel('v'); zlabel('Pattern');
%
%   See also phitheta2uv, uv2phitheta, uv2phithetapat, phitheta2azelpat.

%   Copyright 2011 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,5,nargin);
eml_assert_no_varsize(1:nargin,pat_old,ang1_old,ang2_old,varargin{:});
if nargin > 4
    ang2_new = varargin{2};
    validateattributes(ang2_new,{'double'},{'vector','real','>=',-1,'<=',1},...
        'phitheta2uvpat','V');
else
    ang2_new = -1:0.01:1;
end
if nargin > 3
    ang1_new = varargin{1};
    validateattributes(ang1_new,{'double'},{'vector','real','>=',-1,'<=',1},...
        'phitheta2uvpat','U');
else
    ang1_new = -1:0.01:1;
end

validateattributes(ang1_old,{'double'},{'vector','real','>=',0,'<=',360},...
    'phitheta2uvpat','PHI');
validateattributes(ang2_old,{'double'},{'vector','real','>=',0,'<=',90},...
    'phitheta2uvpat','THETA');
validateattributes(pat_old,{'double'},{'2d','real','size',[numel(ang2_old) numel(ang1_old)]},...
    'phitheta2uvpat','PAT_PHITHETA');

fh = @uv2phitheta;
pat_new = phased.internal.convertpat(fh,pat_old,ang1_old(:),ang2_old(:),ang1_new(:),ang2_new(:),true);



% [EOF]
