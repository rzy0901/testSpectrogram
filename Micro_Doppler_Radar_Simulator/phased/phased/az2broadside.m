function bsd = az2broadside(az,el)
%az2broadside Convert azimuth to broadside angle
%   BSANG = az2broadside(AZ,EL) returns the broadside angle BSANG
%   corresponding to azimuth AZ and elevation EL. Both AZ and EL are in the
%   local coordinate system. The angles are all in degrees. The azimuth
%   must be within [-180 180] and the elevation must be within [-90 90].
%   The returned broadside angle is within [-90 90].
%
%   AZ and EL can be either scalars or vectors. If both of them are
%   vectors, their dimensions must match.
%
%   BSANG = az2broadside(AZ) returns the broadside angle for 0 elevation.
%
%   % Example:
%   %   Calculate the broadside angle corresponding to azimuth 45 degrees
%   %   and elevation 20 degrees.
%
%   bsang = az2broadside(45,20)
%
%   See also phased, broadside2az.

%   Copyright 2010-2011 The MathWorks, Inc.
%     

%#codegen
%#ok<*EMCA>

    phased.internal.narginchk(1,2,nargin);
    if nargin < 2
        el = 0;
    end
    sigdatatypes.validateAngle(az,'az2broadside','AZ',...
        {'vector','>=',-180,'<=',180});
    sigdatatypes.validateAngle(el,'az2broadside','EL',...
        {'vector','>=',-90,'<=',90}); 
    cond = (numel(az)>1) && (numel(el)>1) && any(size(az)~=size(el));
    if cond
        coder.internal.errorIf(cond, ...
             'phased:system:DimensionMismatch','AZ','EL');
    end
    
    
    bsd = asind(sind(az).*cosd(el));


% [EOF]
