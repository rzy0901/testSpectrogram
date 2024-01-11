function SNR = shnidman(Pd,Pfa,N,SwerlingCase)
%shnidman Shnidman's equation
%   SNR = shnidman(Pd,Pfa) returns the required SNR (in dB) of a single
%   received sample to achieve the specified probability of detection, Pd,
%   and probability of false alarm, Pfa, when the square law detector is
%   used. Pd should be between 0.1 and 0.99. The SNR is calculated using
%   Shnidman's equation. By default, the target is assumed to be
%   non-fluctuating.
%
%   SNR = shnidman(Pd,Pfa,N) uses N as the number of samples for
%   noncoherent integration.
%
%   SNR = shnidman(Pd,Pfa,N,SwCase) uses SwCase as the target's Swerling
%   case number. SwCases can be 0, 1, 2, 3, 4, 5. Swerling case 0 or 5
%   refers to non-fluctuating targets.
%
%   In general, the estimated SNR from Shnidman's equation has a looser
%   bound compared to the estimated SNR from Albersheim's equation. For
%   typical values of Pd within [0.1 0.99], Pfa within [1e-9 1e-3], and N
%   within [1 100], the error is less than 0.5 dB.
%
%   % Example:
%   %   Find out the required SNR to achieve probability of detection of
%   %   0.9 and probability of false alarm of 1e-3 using Shnidman's
%   %   equation. Assuming that the target is a Swerling case 1 target and
%   %   10 samples are integrated.
%
%   snr = shnidman(0.9,1e-3,10,1)
%
%   See also phased, albersheim.

%   Copyright 2010-2016 The MathWorks, Inc.

%   Reference
%   [1] Mark Richards, Fundamentals of Radar Signal Processing, Page 337

%#codegen
%#ok<*EMCA>

if nargin < 4
    SwerlingCase = 0;
end

if nargin < 3
    N = 1;
end

localvalidateinput(Pd,Pfa,N,SwerlingCase);

switch SwerlingCase
    case 1
        K = 1;
    case 2
        K = N;
    case 3
        K = 2;
    case 4
        K = 2*N;
    otherwise %case {0,5}
        K = inf;
end

if N < 40
    alpha = 0;
else
    alpha = 1/4;
end

eta = sqrt(-0.8*log(4*Pfa*(1-Pfa)))+sign(Pd-0.5)*sqrt(-0.8*log(4*Pd*(1-Pd)));
Xinf = eta*(eta+2*sqrt(N/2+(alpha-1/4)));
C1 = (((17.7006*Pd-18.4496)*Pd+14.5339)*Pd-3.525)/K;
C2 = 1/K*(exp(27.31*Pd-25.14)+(Pd-0.8)*(0.7*log(1e-5/Pfa)+(2*N-20)/80));

if Pd<=0.872
    Cdb = C1;
else 
    Cdb = C1+C2;
end

SNR = Cdb + pow2db(Xinf/N);

%-------------------------------------------
function localvalidateinput(pd, pfa, N, Sw)

eml_assert_no_varsize(1:nargin, pd, pfa, N, Sw);
validateattributes(pd,{'numeric'},{'positive','nonempty','scalar','>=',.1,'<=',.99},...
    'shnidman','Pd');

validateattributes(pfa,{'numeric'},{'positive','nonempty','scalar','<=',1},...
    'shnidman','Pfa');

validateattributes(N,{'numeric'},{'positive','nonempty','integer','scalar','finite'},...
    'shnidman','N');

validateattributes(Sw,{'numeric'},{'nonempty','scalar','integer','>=',0,'<=',5},...
    'shnidman','SwCase');


% [EOF]
