function epat = privgpuCrossedDipoleAntElemResp(az,el, epat)
%This function is for internal use only. It may be removed in the future.

%PRIVGPUCROSSEDDIPOLEANTELEMRESP - computes the element response of a
%CrossedDipoleAntennaElement for a given az,el (in radians).

%   Copyright 2013 The MathWorks, Inc.

scalefactor = sqrt(3/2);  

%H-response
h = -1i.*cos(az);
%V-response
fh = @(x,y) -cos(y) + 1i.*sin(y).*sin(x);
v = bsxfun(fh, az, el);
epat(:) = scalefactor.*bsxfun(@hypot,h,v);
