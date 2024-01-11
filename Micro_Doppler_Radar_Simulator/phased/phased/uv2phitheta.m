function phitheta = uv2phitheta(uv)
%UV2PHITHETA  Convert angles from u/v format to phi/theta format
%   PhiTheta = uv2azel(UV) converts the u/v value pairs to their
%   corresponding phi/theta angles. UV is a 2-row matrix whose columns are
%   values specified in [u;v] format, where both u and v are between -1 and
%   1. PhiTheta has the same dimension as UV. The columns of PhiTheta are
%   the angles specified in [phi;theta] (in degrees) form.
%
%   Assume the boresight is the X-axis. The phi angle is the angle from the
%   Y-axis to the Z-axis, ranging from 0 and 360 degrees. The theta angle
%   is the angle from X-axis toward the Y-Z plane, ranging from 0 to 90
%   degrees. The value u is defined as sin(theta)*cos(phi) and the value v
%   is defined as sin(theta)*sin(phi).
%
%   % Example: 
%   %   Find the corresponding phi and theta angles for u of 0.5 and v of
%   %   0.
%   
%   phitheta = uv2phitheta([0.5;0])
%
%   See also phitheta2uv, phitheta2uvpat, uv2azel, uv2phithetapat.

%   Copyright 2011-2014 The MathWorks, Inc.

%#codegen 
%#ok<*EMCA>

% u and v between [-1 1], within unit circle
validateattributes(uv,{'double'},{'2d','real','>=',-1,'<=',1},...
    'uv2phitheta','UV');

cond = size(uv,1) == 2;
if ~cond
    coder.internal.assert(cond,'phased:phased:invalidRowNumbers','UV',2);
end

u = uv(1,:);
v = uv(2,:);

cond = all(hypot(u,v)<=1);
if ~cond
    coder.internal.assert(cond,'phased:phased:invalidUV');
end

% phi
phi = atand(v./u);

% take care of quadrants
temp_idx = (u<0);
phi(temp_idx) = phi(temp_idx)+180;
temp_idx = (u>0)&(v<0);
phi(temp_idx) = phi(temp_idx)+360;

% take care of planes
phi((u<0)&(v==0)) = 180;
phi((u==0)&(v<0)) = 270;

% take care of axes
phi((u==0)&(v==0)) = 0;
phi((u==-1)&(v==0)) = 180;
phi((u==0)&(v==-1)) = 270;

% theta
temp = u.^2+v.^2;
temp(abs(temp)>1)=1;
theta = asind(sqrt(temp));

phitheta = [phi;theta];

% [EOF]
