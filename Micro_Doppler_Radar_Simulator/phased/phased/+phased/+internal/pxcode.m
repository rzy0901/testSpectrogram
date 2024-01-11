function x = pxcode(Nchip)
%This function is for internal use only. It may be removed in the future.

%PXCODE Px code
%   X = phased.internal.pxcode(M) returns an M-length Px code in X. M must
%   be a perfect square. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip Px code.
%
%   x = phased.internal.pxcode(16);

%   Copyright 2011 The MathWorks, Inc.

% Assume Nchip is already validated as a perfect square.

%#codegen

Nsqrt = sqrt(Nchip);
chipidx = (0:Nsqrt-1).';
if rem(Nsqrt,2)
    % Nsqrt odd
    phasemat = ((Nsqrt-1)/2-chipidx)*(Nsqrt/2-1-chipidx).';
else
    % Nsqrt even
    temp = (Nsqrt-1)/2-chipidx;
    phasemat = temp*temp.';
end
x = exp(1i*2*pi*phasemat.'/Nsqrt);
x = x(:);



% [EOF]
