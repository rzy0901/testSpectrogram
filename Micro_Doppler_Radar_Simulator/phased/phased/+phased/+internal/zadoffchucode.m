function x = zadoffchucode(Nchip,SeqIdx)
%This function is for internal use only. It may be removed in the future.

%ZADOFFCHUCODE Zadoff-Chu code
%   X = phased.internal.zadoffchucode(M,SEQIDX) returns an M-length
%   Zadoff-Chu code with sequence index SEQIDX in X. M must be a positive
%   integer and SEQIDX must be prime to M. X is a column vector.
%
%   Example:
%   %   Generate a 16-chip Zadoff-Chu code.
%
%   x = phased.internal.zadoffchucode(16,1);

%   Copyright 2011 The MathWorks, Inc.

% Assume it is already validated that Nchip and SeqIdx is prime to each
% other

%#codegen

m = (0:Nchip-1).';
if rem(Nchip,2)
    % Nchip odd
    phasevec = m.*(m+1)/2;
else
    % Nchip even
    phasevec = m.^2/2;
end
phasevec = SeqIdx*phasevec;
x = exp(1i*2*pi*phasevec/Nchip);



% [EOF]
