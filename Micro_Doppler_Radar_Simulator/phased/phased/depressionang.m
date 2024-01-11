function depAng = depressionang(platformHeight,Rclutter,earthmodel,Re)
%depressionang Depression angle of surface target
%   DEPANG = depressionang(H,R) returns the depression angle (in degrees)
%   toward surface targets for a sensor at the height of H (in meters)
%   above the surface. R is the range (in meters) between the surface
%   target and the sensor. Horizontal is considered as 0 degrees depression
%   angle. The computation assumes curved earth model with 4/3 effective
%   earth radius.
%
%   H and R can be either scalars or vectors. If both H and R are vectors,
%   then the dimensions of the two vectors must be the same. R must be
%   greater than or equal to the corresponding H and less than or equal to
%   the horizon range.
%
%   DEPANG = depressionang(H,R,MODEL) specifies the earth model used to
%   compute the depression angle as one of 'Flat' | 'Curved', where the
%   default is 'Curved'. 
%
%   DEPANG = depressionang(H,R,MODEL,RE) specifies the effective earth
%   radius (in meters) as a positive scalar, RE, where the default is 4/3
%   of the earth radius. This input is ignored when you set MODEL to
%   'Flat'.
%
%   % Example:
%   %   Calculate the depression angle for a ground clutter patch that is 
%   %   1000 meters away from the sensor.  Sensor is located on a platform 
%   %   that is 300 meters above the ground.
%
%   depang = depressionang(300,1000)
%
%   See also horizonrange, grazingang.

%   Copyright 2008-2011 The MathWorks, Inc.

%   Reference
%   [1] James Ward, Space-Time Adaptive Processing for Airborne Radar, 1994
%   [2] Long, Radar Reflectivity of Land and Sea, 2001

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(2,4,nargin);

if nargin < 4
    Re = effearthradius;
end

if nargin < 3
    earthmodel = 'Curved';
end
eml_assert_no_varsize(3:4, platformHeight,Rclutter,earthmodel,Re);

sigdatatypes.validateDistance(platformHeight,'depressionang','H',{'vector'});
sigdatatypes.validateDistance(Rclutter,'depressionang','R',{'vector','positive'});
earthmodel = validatestring(earthmodel,{'Flat','Curved'},'depressionang','MODEL');
sigdatatypes.validateDistance(Re,'depressionang','Re',{'scalar','positive'});


cond =  numel(Rclutter)>1 && (numel(platformHeight)>1) && ...
        any(size(Rclutter)~=size(platformHeight));
if cond
    coder.internal.errorIf(cond, ...
         'phased:system:DimensionMismatch','H','R');
end


cond = any(Rclutter(:)<platformHeight(:));
if cond
    coder.internal.errorIf(cond,...
         'phased:phased:expectedGreaterThanOrEqualTo','R','H');
end

cond = earthmodel(1) == 'C' && ...
       any(Rclutter(:)>horizonrange(platformHeight(:),Re));
if cond
    coder.internal.errorIf(cond, ...
         'phased:phased:expectedLessThanHorizonRange','R');
end


if earthmodel(1) == 'C' 
    
    depAng = (platformHeight.^2+2.*platformHeight.*Re+Rclutter.^2)./(2.*Rclutter.*(Re+platformHeight));
    depAng = asind(depAng);    

else
    depAng = asind(platformHeight./Rclutter);
end





% [EOF]
