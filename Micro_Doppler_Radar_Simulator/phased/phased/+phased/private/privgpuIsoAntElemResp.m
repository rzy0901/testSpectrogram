function epat = privgpuIsoAntElemResp(az, el, backbaffled)
%This function is for internal use only. It may be removed in the future.

%PRIVGPUISOANTELEMRESP - computes the element response of a
%IsotropicAntennaElement for a given az,el (in radians).

%   Copyright 2012 The MathWorks, Inc.

if backbaffled,
    if isvector(az),
        az = reshape(az, 1,1,[]);
        el = reshape(el, [], 1,1);
    end
    epat = arrayfun(@gpuIsoResp, az, el);    

else
    
    if isvector(az)
        epat = gpuArray.ones(numel(el), 1, numel(az));
    else
        epat = gpuArray.ones(size(el));
    end
end

end

function resp = gpuIsoResp(az, el)
%Computes the isotropic response
    ninety = pi/2;

    rooteps = sqrt(eps);
    
    tf = (az < (-1*ninety - rooteps)) | ...
         (az > (ninety + rooteps)) & ...
         (abs(abs(el) - ninety) > rooteps);

    if tf
        resp = 0;
    else
        resp = 1;
    end
     
     
end

