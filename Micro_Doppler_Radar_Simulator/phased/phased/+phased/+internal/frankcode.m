function x = frankcode(Nchip)
%This function is for internal use only. It may be removed in the future.

%FRANKCODE Frank code
%   X = phased.internal.frankcode(M) returns an M-length Frank code in X. M
%   must be a perfect square. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip Frank code.
%
%   x = phased.internal.frankcode(16);

%   Copyright 2011 The MathWorks, Inc.

% Assume Nchip is already validated as a perfect square.

%#codegen

Nsqrt = sqrt(Nchip);
chipidx = (0:(Nsqrt-1)).';
phasemat = chipidx*chipidx.';
x = exp(1i*2*pi*phasemat.'/Nsqrt);
x = x(:);



% [EOF]
