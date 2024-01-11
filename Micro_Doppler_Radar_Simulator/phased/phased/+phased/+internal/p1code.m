function x = p1code(Nchip)
%This function is for internal use only. It may be removed in the future.

%P1CODE P1 code
%   X = phased.internal.p1code(M) returns an M-length P1 code in X. M must
%   be a perfect square. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip P1 code.
%
%   x = phased.internal.p1code(16);

%   Copyright 2011 The MathWorks, Inc.

% Assume Nchip is already validated as a perfect square.

%#codegen

Nsqrt = sqrt(Nchip);
chipidx = (0:Nsqrt-1).';
temp = (Nsqrt-1)/2-chipidx;
phasemat = -(repmat(temp.*(chipidx*Nsqrt),1,Nsqrt) + ...
    temp*chipidx.');   
x = exp(1i*2*pi*phasemat.'/Nsqrt);
x = x(:);



% [EOF]
