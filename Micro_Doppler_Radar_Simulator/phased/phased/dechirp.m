function y = dechirp(x,xref)
%dechirp   Dechirp FMCW signal
%   Y = dechirp(X,XREF) returns result of mixing the incoming signal, X,
%   with the reference signal XREF. In an FMCW radar system, X is the
%   received signal and XREF is the transmitted signal.
%
%   X can be an MxN matrix while XREF must be an Mx1 vector. If X is a
%   matrix, then each column in X is an independent received signal and is
%   individually mixed with XREF. The output Y has the same dimensionality
%   as X. Each column in Y contains the mixer output for each corresponding
%   signal in X.
%
%   X and XREF can be baseband complex signals. For a column vector X, the
%   mix operation is defined as
%
%   Y = XREF.*conj(X)
%
%   Note that the order of XREF and X implies that the Doppler shift
%   embedded in X will be negated.
%
%   This function supports single and double precision for input data and
%   arguments. If the input data, X, is single precision, the output data 
%   is single precision. If the input data, X, is double precision, the 
%   output data is double precision. The precision of the output is
%   independent of the precision of the arguments.
%
%   % Example
%   %   Create an FMCW signal, delay it and then dechirp it using the
%   %   original FMCW signal. Plot the spectrum of the dechirped signal.
%
%   hwav = phased.FMCWWaveform;
%   xref = step(hwav);
%   x = [zeros(10,1);xref(1:end-10)];  % delay by 10 samples
%   y = dechirp(x,xref);
%   pwelch(y,hamming(64),[],1024,hwav.SampleRate,'centered');
%
%   See also phased, beat2range, phased.FMCWWaveform, range2beat,
%   rdcoupling, stretchfreq2rng.

%   Copyright 2012 The MathWorks, Inc.

%   References
%   [1] Merrill Skolnik, Introduction to Radar Systems, McGraw-Hill, 1962
%   [2] Philip E. Pace, Detecting and Classifying Low Probability of
%       Intercept Radar, 2nd Ed., Artech House, 2009

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,2,nargin);

validateattributes(x,{'double','single'},...
    {'2d','nonempty','finite','nonnan','nonsparse'},'dechirp','X');
validateattributes(xref,{'double','single'},...
    {'column','nonempty','finite','nonnan','nonsparse'},'dechirp','XREF');

xref = cast(xref,class(x));

cond = size(x,1) ~= size(xref,1);
if cond
    coder.internal.errorIf(cond, ...
         'phased:phased:NumRowsMismatch','X','XREF');
end

y = bsxfun(@times,xref,conj(x));


% [EOF]
