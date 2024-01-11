function r = time2range(t,c)
%time2range     Convert propagation time to propagation distance
%   R = time2range(T) returns the distance propagated, R (in meters), of a
%   signal during the time specified in T (in seconds). The propagation is
%   assumed to be two way, as in a monostatic radar system. T must be
%   positive and R has the same dimensionality as T.
%
%   R = time2range(T,C) specifies the signal propagation speed, C (in
%   meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   % Example:
%   %   Calculate the minimum detectable range for a monostatic radar
%   %   system where the pulse width is set to 2 ms.
%
%   tau = 2e-3; c = 3e8;
%   rmin = time2range(tau,c)
%
%   See also phased, range2time, range2bw.

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

sigdatatypes.validateDuration(t,'time2range','T');
sigdatatypes.validateSpeed(c,'time2range','C',{'positive','scalar'});

r = c*t/2;

% [EOF]
