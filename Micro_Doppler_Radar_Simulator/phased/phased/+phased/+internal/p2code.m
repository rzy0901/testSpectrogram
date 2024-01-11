function x = p2code(Nchip)
%This function is for internal use only. It may be removed in the future.

%P2CODE P2 code
%   X = phased.internal.p2code(M) returns an M-length P2 code in X. M must
%   be an even perfect square. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip P2 code.
%
%   x = phased.internal.p2code(16);

%   Copyright 2011 The MathWorks, Inc.

% Assume Nchip is already validated as a perfect square. 
% Assume Nchip is also validated as a even number

%#codegen
x = phased.internal.pxcode(Nchip);



% [EOF]
