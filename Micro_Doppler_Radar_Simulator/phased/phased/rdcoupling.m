function deltaR = rdcoupling(Fd,beta,c)
%rdcoupling     Range Doppler coupling
%   DR = rdcoupling(FD,SLOPE) returns the range offset, DR (in meters), due
%   to Doppler shift, FD (in Hz), in a linear frequency modulation signal,
%   such as a linear FM pulse or an FMCW signal. DR has the same
%   dimensionality as FD.
%
%   SLOPE specifies the sweep slope (in Hz/s) of the linear frequency
%   modulation as a nonzero scalar.
%
%   R = rdcoupling(FB,SLOPE,C) specifies the signal propagation speed, C
%   (in meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   The range offset is computed based on the following equation:
%
%   DR = -C*FD/(2*SLOPE)
%
%   Note that range offset is defined as the difference between the
%   estimated range and the true range:
%
%   range offset = estimated range - true range
%
%   % Example:
%   %   An FMCW waveform sweeps a band of 3 MHz in 2 ms. The dechirped
%   %   signal of a target return has a beat frequency of 1 kHz. The
%   %   processing of the signal also revealed a Doppler shift of 100 Hz.
%   %   Calculate the true range of the target.
%
%   swslope = 30e6/2e-3;
%   fb = 1e3; fd = 100; c = 3e8;
%   r = beat2range(fb,swslope,c) - rdcoupling(fd,swslope,c)
%
%   See also phased, beat2range, dechirp, phased.FMCWWaveform,
%   phased.LinearFMWaveform, range2beat, stretchfreq2rng.

%   Copyright 2012 The MathWorks, Inc.

%   References
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005
%   [2] David K. Barton, Radar System Analysis and Modeling, Artech House,
%       2005

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,3,nargin);

if nargin < 3
    c = physconst('lightspeed');
end

validateattributes(Fd,{'double'},{'real','nonempty','finite'},...
    'rdcoupling','FD');

validateattributes(beta,{'double'},{'real','finite','nonempty',...
    'scalar','nonzero'},'rdcoupling','SLOPE');

sigdatatypes.validateSpeed(c,'rdcoupling','C',{'positive','scalar'});

deltaR = -c*Fd/(2*beta);

% [EOF]
