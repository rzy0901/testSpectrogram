function P = billingsleyicm(Fd,Fc,wspeed,c)
%billingsleyicm Billingsley's intrinsic clutter motion (ICM) model
%   P = billingsleyicm(FD,FC,WSPEED) calculates the clutter Doppler
%   spectrum shape, P, due to intrinsic clutter motion (ICM) at Doppler
%   frequencies specified in FD (in Hz) using Billingsley's model. FD can
%   be a vector and P has the same dimensionality as FD. FC is the
%   operating frequency (in Hz) of the system and WSPEED is the wind speed
%   (in m/s).
%
%   ICM arises when wind blows on vegetation or other clutter sources. Only
%   a small range of Doppler frequencies contribute significantly to the
%   Doppler spectrum. Typically, these Doppler frequencies correspond to
%   speeds less than 3 m/s.
%
%   P = billingsleyicm(FD,FC,WSPEED,C) specifies the propagation speed (in
%   m/s) in C. The default value is the speed of light.
%
%   % Example:
%   %   Calculate and plot the Doppler spectrum shape predicted by 
%   %   Billingsley's ICM model. Assume the PRF is 2 kHz, the operating
%   %   frequency is 1 GHz and the wind speed is 5 m/s.
%
%   v = -3:0.1:3; fc = 1e9; wspeed = 5; c = 3e8;
%   fd = 2*v/(c/fc);
%   p = billingsleyicm(fd,fc,wspeed);
%   plot(fd,pow2db(p)); 
%   xlabel('Doppler frequency (Hz)'), ylabel('P (dB)');

%   Copyright 2011 The MathWorks, Inc.
%     

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(3,4,nargin);

if nargin < 4
    c = physconst('lightspeed');
end

[c, uratio] = parseInput(Fd,Fc,wspeed,c);    

lambda = c/Fc;
% Billingsley's uses wind speed in status miles per hour and operating
% frequency in GHz
wspeed = wspeed*3600*uratio;
Fc = Fc/1e9;

ZeroFdIdx = find(Fd==0);

alpha = 489.8*(wspeed^-1.55)*(Fc^-1.21);
beta = 1/(0.1048*(log10(wspeed)+0.4147));

P = 1/(1+alpha)*(lambda*beta/4)*exp(-(lambda*beta/2)*abs(Fd));
P(ZeroFdIdx) = P(ZeroFdIdx) + alpha/(1+alpha);


function [c, uratio] = parseInput(Fd,Fc,wspeed,c)

eml_assert_no_varsize(1:nargin, Fd,Fc,wspeed,c);
coder.extrinsic('unitsratio');
uratio = coder.internal.const(unitsratio('sm','m'));

sigdatatypes.validateFrequency(Fc,'billingsleyicm','FC',{'scalar'});
sigdatatypes.validateSpeed(wspeed,'billingsleyicm','WSPEED',{'scalar'});
sigdatatypes.validateSpeed(c,'billingsleyicm','C',{'scalar','positive'});
validateattributes(Fd,{'double'},{'vector','finite','real'},'billingsleyicm',...
    'FD');

% [EOF]
