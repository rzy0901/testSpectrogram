function snr = albersheim(pd, pfa, N)
%albersheim Albersheim's equation
%   SNR = albersheim(PD,PFA) returns the required SNR (in dB) of the single
%   received sample to achieve the specified probability of detection, PD,
%   and probability of false alarm, PFA, when noncoherent detection is
%   used.  The SNR is calculated using Albersheim's equation. Albersheim's
%   equation assumes that the target is non-fluctuating.
%
%   SNR = albersheim(PD,PFA,N) also considers the noncoherent pulse
%   integration of N pulses.  The default value of N is 1.
%
%   The error in the estimated SNR is less than 0.2 dB when PFA is within
%   [1e-7 1e-3], PD is within [0.1 0.9], and N is within [1 8096].
%
%   % Example:
%   %   Determine the required SNR to achieve a probability of detection 
%   %   of 0.9 and probability of false alarm of 1e-3 using Albersheim's
%   %   equation.
%
%   snr = albersheim(0.9,1e-3)
%
%   See also phased, shnidman.

%   Copyright 2007-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, pg 329
%   [2] Merrill Skolnik, Introduction to Radar Systems, pg 49

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(2,3,nargin);

if nargin < 3
    N = 1;
end
eml_assert_no_varsize(1:nargin, pd, pfa, N);
localvalidateinput(pd, pfa, N);

A = log(0.62/pfa);
B = log(pd/(1-pd));

N1 = 1/sqrt(N);
N2 = (6.2+4.54/sqrt(N+0.44))/10;
snr = N1*(A + 0.12*A*B + 1.7*B)^N2;

snr = pow2db(snr);

%-------------------------------------------
function localvalidateinput(pd, pfa, N)

validateattributes(pd,{'numeric'},{'positive','nonempty','scalar','<',1},...
    'albersheim','PD');

validateattributes(pfa,{'numeric'},{'positive','nonempty','scalar','<',1},...
    'albersheim','PFA');

validateattributes(N,{'numeric'},{'positive','nonempty','integer','scalar','finite'},...
    'albersheim','N');


% [EOF]
