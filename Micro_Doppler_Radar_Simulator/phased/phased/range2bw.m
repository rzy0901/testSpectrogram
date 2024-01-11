function bw = range2bw(r,c)
%range2bw   Convert range resolution to required bandwidth
%   BW = range2bw(R) returns the bandwidth, BW (in Hz), needed to
%   distinguish two targets separated by the range specified in R (in
%   meters). Such capability is often referred to as range resolution. The
%   propagation is assumed to be two way, as in a monostatic radar system.
%   R must be positive and BW has the same dimensionality as R.
%
%   BW = range2bw(R,C) specifies the signal propagation speed, C (in
%   meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   % Example:
%   %   Calculate the required pulse width for a rectangular waveform used 
%   %   in a monostatic radar system so that it can achieve a range 
%   %   resolution of 10 m.
%
%   delta_r = 10; c = 3e8;
%   tau = 1/range2bw(delta_r,c)
%
%   See also phased, time2range, range2time.

%   Copyright 2012 The MathWorks, Inc.

%   References
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Edition,
%       McGraw-Hill, 2001

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(1,2,nargin);

if nargin < 2
    c = physconst('lightspeed');
end

sigdatatypes.validateDistance(r,'range2bw','R',{'positive'});
sigdatatypes.validateSpeed(c,'range2bw','C',{'positive','scalar'});

bw = c./(2*r);


% [EOF]
