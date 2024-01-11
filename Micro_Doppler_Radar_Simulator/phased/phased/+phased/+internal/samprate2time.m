function st = samprate2time(sr,N)
%This function is for internal use only. It may be removed in the future.

%SAMPRATE2TIME Translate sample rate to sample time
%   ST = phased.internal.samprate2time(SR,N) converts the signal sample
%   rate (in Hz), SR, to the corresponding sample time (in second), ST, of
%   a length-N signal vector .

%   Copyright 2014 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

st = N/sr;


% [EOF]
