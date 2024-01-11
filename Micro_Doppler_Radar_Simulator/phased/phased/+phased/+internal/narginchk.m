function narginchk(low,high,nargin)

%   Copyright 2009-2011 The MathWorks, Inc.

%#codegen

cond = nargin >= low;
if ~cond
    coder.internal.assert(cond,'MATLAB:narginchk:notEnoughInputs');
end

cond = nargin <= high;
if ~cond
    coder.internal.assert(cond,'MATLAB:narginchk:tooManyInputs');
end

end
