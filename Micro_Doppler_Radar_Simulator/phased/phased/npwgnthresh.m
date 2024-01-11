function T = npwgnthresh(pfa,npulses,dtype,oscale)
%npwgnthresh Detection SNR threshold for white Gaussian noise
%   T = npwgnthresh(PFA) calculates the SNR threshold T (in dB) for
%   detecting a signal with the given probability of false alarm, PFA, in
%   wide sense stationary white Gaussian noise using the Neyman-Pearson
%   (NP) decision rule.
%
%   T = npwgnthresh(PFA,NPULS) specifies NPULS as the number of pulses used
%   in the pulse integration. The default value of NPULS is 1.
%
%   T = npwgnthresh(PFA,NPULS,DTYPE) specifies DTYPE as the type of
%   detection using one of 'coherent' | 'noncoherent' | 'real', where the
%   default is 'noncoherent'. A square law detector is used in noncoherent
%   detection. For a nonfluctuating target, the detection can be performed
%   either coherently or noncoherently. For a fluctuating target, such as
%   Swerling cases 1-4, a noncoherent detection must be performed.
%
%   T = npwgnthresh(PFA,NPULS,DTYPE,OUTSCALE) specifies OUTSCALE as the
%   scale of the output value using one of 'db' | 'linear', where the
%   default is 'db'. When the OUTSCALE is set to 'linear', the returned
%   threshold is represented in amplitude. 
%   
%   Note that when PFA is greater than 0.5 and the DTYPE is set to
%   'coherent' or 'real', the resulting threshold cannot be represented in
%   dB scale.
%
%   % Examples:
%   
%   % Example 1:
%   %   Calculate the SNR threshold that achieves a probability of false
%   %   alarm (Pfa) of 0.01 using a detection type of 'real'.
%
%   snrthreshold = npwgnthresh(0.01, 1,'real');
%
%   % Example 2:
%   %   Verify that the above threshold produces a Pfa of 0.01 by
%   %   constructing 10000 white real Gaussian noise samples and counting
%   %   how many times the sample passes the threshold.
%
%   npower = 1; Ntrial = 10000;
%   noise = sqrt(npower)*randn(1,Ntrial);
%   snrthreshold = npwgnthresh(0.01, 1,'real');
%   threshold = sqrt(npower*db2pow(snrthreshold));
%   Pfa = sum(noise>threshold)/Ntrial
%
%   % Example 3:
%   %   Plot the curve of SNR threshold as a function of probability of
%   %   false alarm (Pfa) for a Swerling case 2 target assuming white
%   %   Gaussian noise. 10 pulses will be noncoherently integrated before
%   %   the detection. Note that for a fluctuating target, one can only use
%   %   noncoherent detection.
%
%   Pfa = 0.1:0.01:0.6;
%   for m = numel(Pfa):-1:1
%       snrthreshold(m) = npwgnthresh(Pfa(m), 10,'noncoherent','linear');
%   end
%   plot(Pfa,snrthreshold); grid on;
%   xlabel('Probability of false alarm'); ylabel('SNR threshold');
%
%   See also phased, rocsnr, rocpfa.

%   Copyright 2007-2012 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing,
%       McGraw-Hill, 2005

%#codegen 
%#ok<*EMCA>

phased.internal.narginchk(1,4,nargin);

if nargin < 4
    oscale = 'db';
end

if nargin < 3
    dtype = 'noncoherent';
end

if nargin < 2
    npulses = 1;
end
eml_assert_no_varsize(1:nargin,pfa,npulses,dtype,oscale);
sigdatatypes.validateProbability(pfa,'npwgnthresh','PFA',...
    {'<',1,'>',0,'scalar'});
sigdatatypes.validateIndex(npulses,'npwgnthresh','NPULS',{'scalar'});
dtype = validatestring(dtype,{'noncoherent','coherent','real'},...
    'npwgnthresh','DTYPE');
oscale = validatestring(oscale,{'db','linear'},...
    'npwgnthresh','OUTSCALE');

cond = ...
    pfa > 0.5 && (dtype(1) == 'c' || dtype(1) == 'r') && (oscale(1) == 'd');
if cond    
    coder.internal.errorIf(cond,'phased:npwgnthresh:invalidScale','PFA',...
        'DTYPE','coherent','real','OUTSCALE','linear');    
end

if dtype(1) == 'c'
    T = erfcinv(2*pfa)*sqrt(npulses);
elseif dtype(1) == 'r' 
    T = erfcinv(2*pfa)*sqrt(2*npulses);
else  
    T = sqrt(gammaincinv(1-pfa,npulses));
end

if oscale(1) == 'd'
    T = mag2db(abs(T));
end


% [EOF]
