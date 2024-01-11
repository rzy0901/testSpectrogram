function epat = privgpuShortDipoleAntElemResp(az,el, axisDir, epat)
%This function is for internal use only. It may be removed in the future.

%PRIVGPUSHORTDIPOLEANTELEMRESP - computes the element response of a
%ShortDipoleAntennaElement for a given az,el (in radians).

%   Copyright 2013 The MathWorks, Inc.

scalefactor = sqrt(3/2);  

%H-response
if strcmp(axisDir,'Y')
    h = -cos(az);
else
    h = epat;
end

%V-response
if strcmp(axisDir,'Y')
    fh = @(x,y) sin(x).*sin(y);
    v = bsxfun(fh, az,el); 
else
    v = -cos(el);
end

epat(:) = scalefactor.*bsxfun(@hypot, h,v);
