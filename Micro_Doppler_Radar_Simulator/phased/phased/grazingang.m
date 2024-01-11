function gAngOut = grazingang(platformHeight,Rclutter,earthmodel,Re)
%grazingang Grazing angle of surface target
%   GRAZANG = grazingang(H,R) returns the grazing angle (in degrees) for a
%   sensor at the height of H (in meters) above the surface. R is the range
%   (in meters) between the surface target and the sensor. The computation
%   assumes curved earth model with 4/3 effective earth radius.
%
%   H and R can be either scalars or vectors. If both H and R are vectors,
%   then the dimensions of the two vectors must be the same. R must be
%   greater than or equal to the corresponding H and less than or equal to
%   the horizon range.
%
%   GRAZANG = grazingang(H,R,MODEL) specifies the earth model used to
%   compute the depression angle as one of 'Flat' | 'Curved', where the
%   default is 'Curved'. 
%
%   GRAZANG = grazingang(H,R,MODEL,RE) specifies the effective earth radius
%   (in meters) as a positive scalar, RE, where the default is 4/3 of the
%   earth radius. This input is ignored when you set Model to 'Flat'.
%
%   % Example:
%   %   Determine the grazing angle of a ground target located 1000 meters 
%   %   away from the sensor.  The sensor is mounted on a platform that is 
%   %   300 meters above the ground.
%
%   gang = grazingang(300,1000)
%
%   See also horizonrange, depressionang.

%   Copyright 2009-2011 The MathWorks, Inc.
%     

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
eml_assert_no_varsize(3:nargin, platformHeight,Rclutter,earthmodel,Re);

sigdatatypes.validateDistance(platformHeight,'grazingang','PlatformHeight',{'vector'});
sigdatatypes.validateDistance(Rclutter,'grazingang','RANGE',{'vector','positive'});
earthmodel = validatestring(earthmodel,{'Flat','Curved'},'grazingang','MODEL');
sigdatatypes.validateDistance(Re,'grazingang','Re',{'scalar','positive'});

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
    gAng = -(Rclutter.^2 + Re^2 - (platformHeight + Re).^2)./(2.*Rclutter.*Re);
    gAng = sign(gAng).*min(abs(gAng),1);
    gAngOut = asind(gAng);
else
    gAngOut = asind(platformHeight./Rclutter);
end

end


% [EOF]
