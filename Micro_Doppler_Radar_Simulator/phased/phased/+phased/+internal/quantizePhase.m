function pq = quantizePhase(p,n,cycle)
%This function is for internal use only. It may be removed in the future.

%quantizePhase  Quantize the phase value
%   PQ = phased.internal.quantizePhase(P,N,PMAX) computes the quantized
%   phase value, PQ, corresponding to the input phase, P. N is the number
%   of quantization bits and PMAX is the maximum phase value represented by
%   the quantization bits. 
%
%   % Example:
%   %   Quantize the 25 degree and 350 degree phase shifts on a 3-bit 
%   %   quantizer. Assume 360 degree is the maximum phase.
%   
%   pq = phased.internal.quantizePhase([25 350],3,360)

%   Copyright 2015 The MathWorks, Inc.

%#codegen

N = 2^n;
pstep = 1/N;
ptemp = p/cycle;
pq = round(mod(ptemp,1)/pstep)*pstep*cycle;
pq(abs(pq-cycle)<eps(cycle)) = 0;
