function x = barkercode(Nchip)
%This function is for internal use only. It may be removed in the future.

%BARKERCODE Barker code
%   X = phased.internal.barkercode(M) returns an M-length Barker code in X.
%   M must be one of the following numbers: 2, 3, 4, 5, 7, 11, and 13. X is
%   a column vector.
%
%   Example:
%   %   Generate a 7-chip Barker code.
%
%   x = phased.internal.barkercode(7);

%   Copyright 2011 The MathWorks, Inc.

% Assume Nchip has been validated to be valid Barker code length

%#codegen

switch Nchip
    case 2
        x = [1 -1];
    case 3
        x = [1 1 -1];
    case 4
        x = [1 1 -1 1];
    case 5
        x = [1 1 1 -1 1];
    case 7
        x = [1 1 1 -1 -1 1 -1];
    case 11
        x = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
    case 13
        x = [1 1 1 1 1 -1 -1 1 1 -1 1 -1 1];
end

x = x(:);


% [EOF]
