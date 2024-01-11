function epat = privgpuCosAntElemResp(az, el, azPow, elPow)
%This function is for internal use only. It may be removed in the future.

%PRIVGPUCOSANTELEMRESP - computes the element response of a
%CosineAntennaElement for a given az,el (in radians) and Azimuth Power
%(azPow) and Elevation Power (elPow).

%   Copyright 2012 The MathWorks, Inc.


if isvector(az),
    az = reshape(az, 1,1,[]);
    el = reshape(el, [], 1,1);
end

epat = arrayfun(@gpuCosResp, az, el, azPow, elPow);

end

function resp = gpuCosResp(az, el, azPow, elPow)
%Computes the cosine response
    ninety = pi/2;

    if (az >= -1*ninety) && (az <= ninety)
        resp = (cos(az).^azPow) .* (cos(el).^elPow);
    else
        resp = 0;
    end
end

