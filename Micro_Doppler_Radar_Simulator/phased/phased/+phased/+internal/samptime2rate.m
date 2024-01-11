function sr = samptime2rate(st,N)
%This function is for internal use only. It may be removed in the future.

%SAMPTIME2RATE Translate sample time to sample rate
%   SR = phased.internal.samptime2rate(ST,N) converts the sample time (in
%   second), ST, of a length-N signal vector to the corresponding sample
%   rate (in Hz), SR, for each signal sample.

%   Copyright 2014 The MathWorks, Inc.

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

sr = N/st;


% [EOF]
