function beatfreq = range2beat(rng,beta,c)
%range2beat     Convert range to beat frequency
%   FB = range2beat(R,SLOPE) converts the range, R (in meters), of a
%   dechirped linear FMCW signal to its corresponding range, beat
%   frequency, FB (in Hz). FB has the same dimensionality as R.
%
%   SLOPE specifies the sweep slope (in Hz/sec) of the FMCW signal as a
%   nonzero scalar. 
%
%   FB = range2beat(R,SLOPE,C) specifies the signal propagation speed, C
%   (in meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   The beat frequency is computed based on the following equation:
%
%   FB = (2*SLOPE*R)/C
%
%   % Example:
%   %   An up sweep FMCW waveform is designed to detect a target as far as
%   %   18 km. The waveform sweeps a 300 MHz band within 1 ms. Calculate
%   %   the maximum beat frequency in the received signal assuming the
%   %   target is not in motion.
%
%   swslope = 300e6/1e-3;
%   rmax = 18e3; c = 3e8;
%   fb = range2beat(rmax,swslope,c)
%
%   See also phased, beat2range, dechirp, phased.FMCWWaveform, rdcoupling,
%   stretchfreq2rng.

%   Copyright 2012 The MathWorks, Inc.

%   References
%   [1] Merrill Skolnik, Introduction to Radar Systems, McGraw-Hill, 1962
%   [2] Philip E. Pace, Detecting and Classifying Low Probability of
%       Intercept Radar, 2nd Ed., Artech House, 2009

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,3,nargin);

if nargin < 3
    c = physconst('lightspeed');
end

sigdatatypes.validateDistance(rng,'range2beat','R');

validateattributes(beta,{'double'},{'real','finite','nonempty',...
    'scalar','nonzero'},'range2beat','SLOPE');

sigdatatypes.validateSpeed(c,'range2beat','C',{'positive','scalar'});


beatfreq = 2*rng/c*beta;

% [EOF]
