function rng = beat2range(beatfreq,beta,c)
%beat2range     Convert beat frequency to range
%   R = beat2range(FB,SLOPE) converts the beat frequency, FB (in Hz), of a
%   dechirped linear FMCW signal to its corresponding range, R (in
%   meters).
%
%   FB can be either an Mx1 vector or an Mx2 matrix. If FB is a vector, the
%   FMCW signal performs either an up sweep or a down sweep. Each entry in
%   FB is a beat frequency in the dechirped signal. If FB is a matrix, the
%   FMCW signal performs a triangle sweep. Each row of FB represents a pair
%   of beat frequencies in the dechirped signal in the form of [UpSweepBeat
%   DownSweepBeat]. R is an Mx1 vector whose rows are the ranges
%   corresponding to the beat frequencies specified in the rows of FB.
%
%   SLOPE specifies the sweep slope (in Hz/s) of the FMCW signal as a
%   nonzero scalar. If the FMCW performs a triangle sweep, SLOPE is the
%   sweep slope of the up sweep half.
%
%   R = beat2range(FB,SLOPE,C) specifies the signal propagation speed, C
%   (in meters/second), as a positive scalar. The default value of C is the
%   speed of light.
%
%   The range is computed based on the following equation:
%
%   R = C*FB/(2*SLOPE)
%
%   If FMCW signal employs triangle sweep, the beat frequency corresponding
%   to the range is derived from each pair of up sweep/down sweep beat
%   frequencies using (UpSweepBeat-DownSweepBeat)/2.
%
%   This function supports single and double precision for input data and
%   arguments. If the input data, FB, is single precision, the output data 
%   is single precision. If the input data, FB, is double precision, the 
%   output data is double precision. The precision of the output is
%   independent of the precision of the arguments.
%
%   % Example:
%   %   An FMCW waveform sweeps a band of 3 MHz in 2 ms. The dechirped
%   %   signal of a target return has a beat frequency of 1 kHz. Find the
%   %   range of the target.
%
%   swslope = 3e6/2e-3;
%   fb = 1e3; c = 3e8;
%   r = beat2range(fb,swslope,c)
%
%   See also phased, dechirp, phased.FMCWWaveform, range2beat, rdcoupling,
%   stretchfreq2rng.

%   Copyright 2012-2015 The MathWorks, Inc.

%   References
%   [1] Merrill Skolnik, Introduction to Radar Systems, McGraw-Hill, 1962
%   [2] Philip E. Pace, Detecting and Classifying Low Probability of
%       Intercept Radar, 2nd Ed., Artech House, 2009

%#codegen
%#ok<*EMCA>

classtouse=class(beatfreq);
phased.internal.narginchk(2,3,nargin);

if nargin < 3
    c = physconst('lightspeed');
end
eml_assert_no_varsize(2:3, beatfreq,beta,c);

validateattributes(beatfreq,{'double','single'},{'real','nonempty','finite','2d'},...
    'beat2range','FB');
ncol_fb = size(beatfreq,2);
cond = ncol_fb > 2;
if cond
    coder.internal.errorIf(cond, ...
         'phased:phased:tooManyColumns','FB',2);
end

if ncol_fb == 1
    validateattributes(beta,{'double','single'},{'real','finite','nonempty',...
        'scalar','nonzero'},'beat2range','SLOPE');
else
    validateattributes(beta,{'double','single'},{'real','finite','nonempty',...
        'scalar','positive'},'beat2range','SLOPE');
end

sigdatatypes.validateSpeed(c,'beat2range','C',{'double','single'},{'positive','scalar'});

c=cast(c,classtouse);
beta=cast(beta,classtouse);

if ncol_fb == 1
    beatfreq_rng = beatfreq;
else
    beatfreq_rng = beatfreq*[1;-1]/2;
end

rng = c*beatfreq_rng/beta/2;

% [EOF]
