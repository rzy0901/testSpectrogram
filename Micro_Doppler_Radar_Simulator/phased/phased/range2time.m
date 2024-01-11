function t = range2time(r,c)
%range2time     Convert propagation distance to propagation time
%   T = range2time(R) returns the time, T (in seconds), that a signal takes
%   to propagate given ranges specified in R (in meters). The propagation
%   is assumed to be two way, as in a monostatic radar system. R must be
%   nonnegative and T has the same dimensionality as R.
%
%   T = range2time(R,C) specifies the signal propagation speed, C (in
%   meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   % Example:
%   %   Calculate the required PRF for a monostatic radar system so that it
%   %   can have a maximum unambiguous range of 15 km.
%
%   rmax = 15e3; c = 3e8;
%   prf = 1/range2time(rmax,c)
%
%   See also phased, time2range, range2bw.

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

sigdatatypes.validateDistance(r,'range2time','R');
sigdatatypes.validateSpeed(c,'range2time','C',{'positive','scalar'});

t = 2*r/c;


% [EOF]
