function uv = phitheta2uv(phitheta)
%PHITHETA2UV  Convert angles from phi/theta format to u/v format
%   UV = phitheta2uv(PhiTheta) converts the phi/theta angle pairs to their
%   corresponding u/v values. PhiTheta is a 2-row matrix whose columns are
%   angles specified in [phi;theta] (in degrees) form. UV has the same
%   dimension as PhiTheta. The columns of UV are angles specified in [u;v]
%   format.
%
%   Assume the boresight is the X-axis. The phi angle is the angle from the
%   Y-axis to the Z-axis, ranging from 0 and 360 degrees. The theta angle
%   is the angle from X-axis toward the Y-Z plane, ranging from 0 to 90
%   degrees. The value u is defined as sin(theta)*cos(phi) and the value v
%   is defined as sin(theta)*sin(phi).
%
%   % Example: 
%   %   Find the corresponding u and v representation for 0 degrees phi and
%   %   30 degrees theta.
%   
%   uv = phitheta2uv([0;30])
%
%   See also phitheta2uvpat, phitheta2azel, uv2phitheta, uv2phithetapat.

%   Copyright 2011-2014 The MathWorks, Inc.

%#codegen 
%#ok<*EMCA>

validateattributes(phitheta,{'double'},{'2d','real','nonnegative'},...
    'phitheta2uv','PhiTheta');

cond = size(phitheta,1)==2;
if ~cond
    coder.internal.assert(cond,'phased:phased:invalidRowNumbers','PhiTheta',2);
end

phi = phitheta(1,:);
theta = phitheta(2,:);

% theta should be limited between [0 90]
sigdatatypes.validateAngle(phi,'phitheta2uv','angle phi',...
    {'>=',0,'<=',360});
sigdatatypes.validateAngle(theta,'phitheta2uv','angle theta',...
    {'>=',0,'<=',90});

u = sind(theta).*cosd(phi);
v = sind(theta).*sind(phi);

uv = [u;v];


% [EOF]
