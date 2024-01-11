function [pat_new,ang1_new,ang2_new] = uv2phithetapat(pat_old,ang1_old,ang2_old,varargin)
%UV2PHITHETAPAT Convert radiation pattern from u/v to phi/theta format
%   PAT_PHITHETA = uv2phithetapat(PAT_UV,U,V) converts the radiation
%   pattern PAT_UV represented in u/v space to phi/theta space. U and V are
%   length-P and length-Q vectors, respectively, representing the u and v
%   values at which the original pattern PAT_UV is sampled. Both U and V
%   must be between -1 and 1. PAT_UV is a QxP matrix.
%
%   PAT_PHITHETA is the converted pattern in phi/theta space, covering phi
%   values from 0 to 360 and theta values from 0 to 90 (in degrees). The
%   pattern is uniformly sampled with a step size of 1 in both phi and
%   theta.
%
%   PAT_PHITHETA = uv2phithetapat(...,PHI,THETA) specifies the grid at
%   which PAT_PHITHETA is sampled in phi/theta space (in degrees). PHI
%   and THETA are length-L and length-M vectors, respectively. PHI must be
%   between 0 and 360 and THETA must be between 0 and 90. PAT_PHITHETA
%   is an MxL matrix.
%
%   [PAT_PHITHETA,PHI,THETA] = uv2phithetapat(...) also returns the grid
%   PHI and THETA in phi/theta space.
%
%   % Example:
%   %   Convert a radiation pattern in u/v space to phi/theta space.
%
%   u = -1:0.01:1; v = -1:0.01:1; pat_uv = ones(201,201);
%   [pat_phitheta,phi,theta] = uv2phithetapat(pat_uv,u,v);
%   surf(phi,theta,pat_phitheta); 
%   xlabel('phi (degrees)'); ylabel('theta (degrees)'); zlabel('Pattern');
%
%   See also uv2phitheta, phitheta2uv, phitheta2uvpat, uv2azelpat.

%   Copyright 2011 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

phased.internal.narginchk(3,5,nargin);
eml_assert_no_varsize(1:nargin,pat_old,ang1_old,ang2_old,varargin{:});
if nargin > 4
    ang2_new = varargin{2};
    validateattributes(ang2_new,{'double'},{'vector','real','>=',0,'<=',90},...
        'uv2phithetapat','THETA');
else
    ang2_new = 0:90;
end
if nargin > 3
    ang1_new = varargin{1};
    validateattributes(ang1_new,{'double'},{'vector','real','>=',0,'<=',360},...
        'uv2phithetapat','PHI');
else
    ang1_new = 0:360;
end

validateattributes(ang1_old,{'double'},{'vector','real','>=',-1,'<=',1},...
    'uv2phithetapat','U');
validateattributes(ang2_old,{'double'},{'vector','real','>=',-1,'<=',1},...
    'uv2phithetapat','V');
validateattributes(pat_old,{'double'},{'2d','real','size',[numel(ang2_old) numel(ang1_old)]},...
    'uv2phithetapat','PAT_UV');

fh = @phitheta2uv;
pat_new = phased.internal.convertpat(fh,pat_old,ang1_old(:),ang2_old(:),ang1_new(:),ang2_new(:),false);


% [EOF]
