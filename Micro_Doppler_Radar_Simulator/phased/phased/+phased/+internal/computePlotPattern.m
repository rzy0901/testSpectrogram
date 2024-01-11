function response = computePlotPattern(response,normalizeflag,unit)
%This function is for internal use only. It may be removed in the future.

%computePlotPattern Compute the pattern to be plot
%   RESP_PLOT = computePlotPattern(RESP,NFLAG,UNIT) returns the response
%   pattern to be plot, RESP_PLOT based on the magnitude pattern, RESP, and
%   configurations such as normalization flag, NFLAG, and the pattern plot
%   unit, UNIT.

%   Copyright 2014 The MathWorks, Inc.

    if normalizeflag
        % response is already magnitude, so no abs needed
        maxresponse = max(max(response));
        if maxresponse ~= 0
            response = response/maxresponse;
        end
    end

    switch unit
        case 'mag',  % already mag.

        case {'pow','power'},
            response = response.^2;

        case 'db',
            response = db(response);
    end

end
