function r = stretchfreq2rng(f,beta,r0,c)
%stretchfreq2rng Convert frequency offset from stretch processing to range 
%   R = stretchfreq2rng(FREQ,SLOPE,REFRNG) returns the range (in meters),
%   R, corresponding to the frequency offset (in Hz), FREQ. FREQ is
%   obtained through stretch processing with a reference range (in meters)
%   of REFRNG and a sweeping slope (in Hz/s) of SLOPE. FREQ can be a vector
%   while both REFRNG and SLOPE are scalars. R has the same dimensions as
%   FREQ.
%
%   R = stretchfreq2rng(...,V) specifies the propagation speed (in m/s), V,
%   as a positive scalar.
%
%   % Example:
%   %   Calculate the range corresponding to a frequency offset of 2 kHz
%   %   obtained from stretch processing. Assume the reference range is
%   %   5000 meters and the sweeping slope is 2 GHz/s.
%
%   r = stretchfreq2rng(2e3,2e9,5000)
%
%   See also phased, ambgfun, beat2range, phased.LinearFMWaveform,
%   phased.StretchProcessor, range2beat, rdcoupling.

%   Copyright 2011-2012 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(3,4,nargin);

if nargin < 4
    c = physconst('lightspeed');
end

validateattributes(f,{'double'},{'vector','real','finite'},...
    'stretchfreq2rng','FREQ');
validateattributes(beta,{'double'},{'scalar','real','finite','nonzero'},...
    'stretchfreq2rng','SLOPE');
sigdatatypes.validateDistance(r0,'stretchfreq2rng','REFRNG',...
    {'scalar','positive'});
sigdatatypes.validateSpeed(c,'stretchfreq2rng','V',{'scalar','positive'});

delta_r = -c*f/(2*beta);
r = r0+delta_r;

% [EOF]
