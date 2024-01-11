function x = p3code(Nchip)
%This function is for internal use only. It may be removed in the future.

%P3CODE P3 code
%   X = phased.internal.p3code(M) returns an M-length P3 code in X. M must
%   be a positive integer. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip P3 code.
%
%   x = phased.internal.p3code(16);

%   Copyright 2011 The MathWorks, Inc.

%#codegen

m = (0:Nchip-1).';
phasevec = m.^2/2;
x = exp(1i*2*pi*phasevec/Nchip);



% [EOF]
