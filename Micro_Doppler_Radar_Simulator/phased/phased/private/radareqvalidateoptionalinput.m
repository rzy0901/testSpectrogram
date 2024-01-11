function radareqvalidateoptionalinput(funName,RCS,Ts,Gain,Loss)
%radareqvalidateoptionalinput Validate optional inputs

%   Copyright 2007-2010 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

validateattributes(RCS,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'RCS');

validateattributes(Ts,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'Ts');

validateattributes(Gain,{'double'},{'real','nonempty'},...
    funName,'Gain');

cond = numel(Gain) <= 2;
if ~cond
    coder.internal.assert(cond,'phased:radareq:invalidDimension','Gain');
end

validateattributes(Loss,{'double'},{'real','nonempty','scalar'},...
    funName,'Loss');

% [EOF]
