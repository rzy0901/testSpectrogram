function y = cwgn(N0,M,N,rs)
%This function is for internal use only. It may be removed in the future.

%CWGN     Complex white Gaussian noise
%   y = phased.internal.cwgn(Pn) returns a complex white Gaussian noise
%   sample, y, whose noise power is specified in N0 (in Watts).
%
%   y = phased.internal.cwgn(Pn, M) returns the noise samples in an Mx1
%   vector where M must be a positive integer.
%
%   y = phased.internal.cwgn(Pn, M, N) returns the noise samples in an
%   M-by-N array.
%
%   y = phased.internal.cwgn(..., RS) uses the random stream specified in
%   RS.
%
%   Example:
%   % Create 10 samples of complex white Gaussian noise whose noise power
%   % is 2 watts.
%   x = phased.internal.cwgn(2,10)

%   Copyright 2010 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

if nargin < 3
    N = 1;
end

if nargin < 2
    M = 1;
end

sigdatatypes.validateIndex(M,'phased.internal.cwgn','M',{'scalar'});
sigdatatypes.validateIndex(N,'phased.internal.cwgn','N',{'scalar'});
sigdatatypes.validatePower(N0,'phased.internal.cwgn','Pn',{'scalar'});

UseGlobalStream = true;
if nargin > 3
    validateattributes(rs,{'RandStream'},{'scalar'},...
        'phased.internal.cwgn','RS')
    UseGlobalStream = false;
end

if UseGlobalStream
    y = sqrt(N0/2)*complex(randn(M,N),randn(M,N));
else
    y = sqrt(N0/2)*complex(randn(rs,M,N),randn(rs,M,N));
end

% [EOF]
