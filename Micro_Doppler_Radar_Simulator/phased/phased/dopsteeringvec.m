function dpstv = dopsteeringvec(dp,numpulses,prf)
%dopsteeringvec Steering vector for Doppler
%   DSTV = dopsteeringvec(DOP,NUMPULSES) returns the steering vector in the
%   time domain corresponding to the Doppler frequency DOP (in Hz). The
%   number of pulses is specified in NUMPULSES and the pulse repetition
%   frequency (PRF) is assumed to be 1 Hz. If DOP is a row vector, DSTV is
%   a matrix with each column containing the steering vector for the
%   corresponding element of DOP.
%
%   DSTV = dopsteeringvec(...,PRF) specifies the pulse repetition frequency
%   PRF (in Hz).
%
%   % Example:
%   %   Calculate the steering vector corresponding to a Doppler frequency
%   %   of 200 Hz, assuming there are 10 pulses and the PRF is 1 kHz.
%
%   dstv = dopsteeringvec(200,10,1000)
%
%   See also phased, speed2dop, dop2speed.

%   Copyright 2008-2010 The MathWorks, Inc.

%   Reference
%   [1] J. R. Guerci, Space-Time Adaptive Processing for Radar, Artech
%       House, 2003

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(2,3,nargin);
if nargin < 3
    prf = 1;
end
eml_assert_no_varsize(2:3, dp,numpulses,prf);
sigdatatypes.validateIndex(numpulses,'dopsteeringvec','NUMPULSES',{'scalar'});
sigdatatypes.validateFrequency(prf,'dopsteeringvec','PRF',{'scalar'});
validateattributes(dp, {'numeric'},...
    {'finite','nonnan','real','vector','>=',-prf/2,'<=',prf/2},...
    'dopsteeringvec','DOP');

% dp could be a vector, make it a row
dp = dp(:).';  
dpstv = exp(1i*2*pi*(0:numpulses-1)'*dp/prf);



% [EOF]
