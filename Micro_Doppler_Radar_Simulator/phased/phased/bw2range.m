function r = bw2range(bw,varargin)
%bw2range   Convert bandwidth to range resolution
%   R = bw2range(BW) returns the range resolution, R (in meters), for the
%   bandwidth specified in BW (in Hz).  The propagation is assumed to be
%   two way, as in a monostatic radar system. BW must be positive and R has
%   the same dimensionality as BW.
%
%   R = bw2range(BW,C) specifies the signal propagation speed, C (in
%   meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   % Example:
%   %   Calculate the range resolution of a monostatic radar operating with
%   %   a system bandwidth of 150 MHz.
%
%   bw = 150e6; c = 3e8;
%   res = bw2range(bw,c)
%
%   See also phased, range2bw, time2range, range2time.

%   Copyright 2016 The MathWorks, Inc.

%   References
%   [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Edition,
%       McGraw-Hill, 2001

%#codegen
%#ok<*EMCA>

narginchk(1,2);

if nargin < 2
    c = physconst('lightspeed');
else
    c = varargin{1};
end

sigdatatypes.validateDistance(bw,'bw2range','BW',{'positive'});
sigdatatypes.validateSpeed(c,'bw2range','C',{'positive','scalar'});

r = c./(2*bw);
