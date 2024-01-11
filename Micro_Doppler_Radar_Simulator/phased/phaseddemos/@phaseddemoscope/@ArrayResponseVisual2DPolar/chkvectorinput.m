 function [angles, response] = chkvectorinput(this)

%   Copyright 2011 The MathWorks, Inc.

    Nresp = validateinputs(this);
    % Loop through all the responses.
    for cnt = 1:Nresp,
        angles(:,cnt) = this(cnt).Angle; %#ok<AGROW>  % the property angle is not included here
        response(:,cnt) = this(cnt).Pattern;  %#ok<AGROW>
    end
 end
    
function [Nresp] = validateinputs(this)
%VALIDATEINPUTS Validate inputs.
%   Validates inputs, and pre-allocates memory.

Nresp = numel(this);
dataLen = numel(this(1).Pattern); % Cache length of 1 data set to compare against.

% Validate inputs and error if necessary.
validateresponses(this,Nresp,dataLen);
end

function validateresponses(this,Nresp,dataLen)
%VALIDATERESPONSES Validate the array responses.
%   Verify that all the responses are the same length.
for k = 2:Nresp,
    validateattributes(this(k).Pattern,{'double'},{'numel',dataLen},...
        '','array response');
end
end