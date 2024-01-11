function r = normalizeIntegratedPower(resp,intpowresp_in,powflag)
%This function is for internal use only. It may be removed in the future.

%normalizeIntegratedPower Compute the field/power normalized to total
%radiated power
%   R = normalizeIntegratedPower(RESP,INTRESP,ISPOWER) returns the response
%   pattern, RESP, normalized by the total power, INTRESP.
%
%   If ISPOWER is true, RESP is a power pattern, the normalization is given
%   by
%
%       R = 4*pi*RESP/INTRESP
%
%   If ISPOWER is false, RESP is a field pattern, the normalization is
%   given by
%
%       R = 2*sqrt(pi)*RESP/sqrt(INTRESP)
%
%   This calculation is useful in computing directivities. It handles the
%   case where the field/power pattern are all zeros, e.g., the Horizontal
%   pattern for a vertical dipole.

%   Copyright 2014 The MathWorks, Inc.

%#codegen

if powflag
    if numel(size(resp)) == 3
        if isscalar(intpowresp_in)
            intpowresp = repmat(intpowresp_in,1,size(resp,3));
        else
            intpowresp = intpowresp_in;
        end
        r = bsxfun(@rdivide,resp,permute(intpowresp,[1 3 2]));
        r(:,:,intpowresp==0) = 0;
    else
        if isscalar(intpowresp_in)
            intpowresp = repmat(intpowresp_in,1,size(resp,2));
        else
            intpowresp = intpowresp_in;
        end
        r = bsxfun(@rdivide,resp,intpowresp);
        r(:,intpowresp==0) = 0;
    end
    r = 4*pi*r;
else
    if numel(size(resp)) == 3
        if isscalar(intpowresp_in)
            intpowresp = repmat(intpowresp_in,1,size(resp,3));
        else
            intpowresp = intpowresp_in;
        end
        r = bsxfun(@rdivide,resp,permute(sqrt(intpowresp),[1 3 2]));
        r(:,:,intpowresp==0) = 0;
    else
        if isscalar(intpowresp_in)
            intpowresp = repmat(intpowresp_in,1,size(resp,2));
        else
            intpowresp = intpowresp_in;
        end
        r = bsxfun(@rdivide,resp,sqrt(intpowresp));
        r(:,intpowresp==0) = 0;
    end
    r = 2*sqrt(pi)*r;
end
    
    