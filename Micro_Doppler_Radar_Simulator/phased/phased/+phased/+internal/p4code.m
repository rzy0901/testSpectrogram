function x = p4code(Nchip)
%This function is for internal use only. It may be removed in the future.

%P4CODE P4 code
%   X = phased.internal.p4code(M) returns an M-length P4 code in X. M must
%   be a positive integer. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip P4 code.
%
%   x = phased.internal.p4code(16);

%   Copyright 2011 The MathWorks, Inc.

%#codegen

m = (0:Nchip-1).';
phasevec = m.*(m-Nchip)/2;
x = exp(1i*2*pi*phasevec/Nchip);



% [EOF]
