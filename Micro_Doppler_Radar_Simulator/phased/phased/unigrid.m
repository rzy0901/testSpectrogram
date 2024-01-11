function y = unigrid(a,step,b,option)
%unigrid Generate a uniform grid
%   Y = unigrid(t1,STEP,t2) returns a uniformly sampled grid Y from closed
%   interval [t1,t2], starting from t1. STEP specifies the step size. This
%   is the same as calling t1:STEP:t2.
%
%   Y = unigrid(...,OPTION) specifies the interval openness. OPTION can be:
%      '[]':  Close at both start and end. (Default)
%      '[)':  Close at start and open at end.
%
%   Note that closed end does not necessarily mean Y will have a sample on
%   t2. Whether t2 is present in Y also depends on the value of STEP.
%
%   % Example:
%   %   Generate a grid on [1:2) with step size 0.1.
%
%   y = unigrid(1,0.1,2,'[)')
%
%   See also phased, val2ind.

%   Copyright 2010-2018 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

    phased.internal.narginchk(3,4,nargin);

    if nargin < 4
        option = '[]';
    end

    validateattributes(a,{'numeric'},{'nonempty','scalar','finite'},...
        'unigrid','t1');
    validateattributes(step,{'numeric'},{'nonempty','scalar','finite'},...
        'unigrid','STEP');
    validateattributes(b,{'numeric'},{'nonempty','scalar','finite'},...
        'unigrid','t2');
    option = validatestring(option,{'[]','[)'},'unigrid','OPTION');
    tmp = a:step:b;
    if ~isempty(tmp) && tmp(end) == b && option(2) ==  ')'
        y = tmp(1:end-1);
    else
        y = tmp;
    end


