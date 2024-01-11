function az = broadside2az(bsd,el)
%broadside2az Convert broadside angle to azimuth
%   AZ = broadside2az(BSANG,EL) returns the azimuth AZ corresponding to the
%   broadside angle BSANG and elevation EL. Both AZ and EL are in the local
%   coordinate system. The angles are all in degrees. The broadside angle
%   and elevation must be within [-90 90]. The returned azimuth is always
%   within [-90 90].
%
%   BSANG and EL can be either scalars or vectors. If both of them are
%   vectors, their dimensions must match.
%
%   AZ = broadside2az(BSANG) returns the azimuth for 0 elevation.
%
%   This function supports single and double precision for inputs BSANG
%   and EL. If any of the inputs is single precision, then the output is
%   single precision.
%
%   % Example:
%   %   Calculate the azimuth corresponding to broadside angle 45 degrees
%   %   and elevation 20 degrees.
%
%   az = broadside2az(45,20)
%
%   See also phased, az2broadside.

%   Copyright 2010-2018 The MathWorks, Inc.
%     

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,2,nargin);
if nargin < 2
    el = 0;
end

sigdatatypes.validateAngle(bsd,'broadside2az','BSANG',...
                           {'double','single'},{'vector','>=',-90,'<=',90});
sigdatatypes.validateAngle(el,'broadside2az','EL',...
                           {'double','single'},{'vector','>=',-90,'<=',90});
cond = (numel(bsd)>1) ...
                           && (numel(el)>1) ...
                           && any(size(bsd)~=size(el));
if cond
    coder.internal.errorIf(cond, ...
                           'phased:system:DimensionMismatch','BSANG','EL');
end

cond = any(abs(el(:))>90-abs(bsd(:)));
if cond
    coder.internal.errorIf(cond, ...
                           'phased:broadside2az:InvalidElevation');
end
    
    az = asind(min(max(sind(bsd)./cosd(el),-1),1));

% [EOF]
