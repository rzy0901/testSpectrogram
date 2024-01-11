function freq = subbandCenterFrequency(fc,fs,Nbands)
%This function is for internal use only. It may be removed in the future.

%subbandCenterFrequency     Compute center frequencies for subbands
%   FREQ = phased.internal.subbandCenterFrequency(FC,FS,NB) returns the
%   center frequencies FREQ as a column vector. FC is the carrier
%   frequency, FS is the sample rate and NB is the number of bands. The
%   function divides [FC-FS/2,FC+FS/2] into NB equal-width subbands.
%
%   Note FREQ is returned corresponding to the result of fft. If one wants
%   to order FREQ from smallest to the greatest, one may want to use
%   fftshift on FREQ.
%
%   % Example
%   %   Compute center frequencies for subbands
%   
%   fc = 3e8; fs = 3e6; nb = 64;
%   freq = phased.internal.subbandCenterFrequency(fc,fs,nb)
%
%   See also phased, phased.internal.SubbandDivider,
%   phased.internal.SubbandCombiner.

%   Copyright 2015 The MathWorks, Inc.

%#ok<*EMCA>
%#codegen

binBound = floor(Nbands/2);
freqBinsR = ((1:Nbands)-binBound-1)*fs/Nbands;
freqBinsC = freqBinsR(:);
freqBins = ifftshift(freqBinsC);
freq = freqBins+fc;

end