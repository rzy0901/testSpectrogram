function SNRdb = radareqsnr(lambda,RNG,Pt,tau,varargin)
%radareqsnr Radar equation to estimate SNR
%   SNR = radareqsnr(LAMBDA,RANGE,Pt,Tau) estimates the signal-to-noise-
%   ratio, SNR (dB), based on the wavelength, LAMBDA (meters), range, RANGE
%   (m), peak power, Pt (Watts), and the pulse width, Tau (seconds).
%
%   For a bistatic radar, specify the range as a two-element row vector
%   representing the range from the transmitter to the target and from the
%   target to the receiver, RANGE = [TxRng RxRng].
%
%   radareqsnr(...,'RCS',RCS) specifies the target's radar cross section
%   RCS in square meters. If omitted RCS = 1 m^2. The function assumes a
%   non-fluctuating target (Swerling case 0).
%
%   radareqsnr(...,'Ts',Ts) specifies the system noise temperature Ts
%   in kelvin (K). If omitted Ts = 290 K.
%
%   radareqsnr(...,'Gain',G) specifies the antenna gain G in decibels (dB).
%   If omitted G = 20 dB. For a bistatic radar, specify the gain as a two
%   element row vector to represent the transmit antenna and receive
%   antenna gain, G = [TxGain RxGain].
%
%   radareqsnr(...,'Loss',L) specifies system losses L in decibels (dB). If
%   omitted L = 0 dB.
%
%   % Example:
%   %   Estimate the SNR required for a radar operating at a frequency of
%   %   5.6 GHz and peak power Pt of 1.5 MW to detect a target at 86.2 km.
%   %   Assume an RCS of 0.1 m^2 and rectangular waveform with a pulse 
%   %   width of 0.2 microseconds. The antenna gain is 45 dB. Assume no
%   %   losses.
%
%   lambda = 3e8/5.6e9; RNG = 85800; Pt = 1.5e6; 
%   tau = 0.2e-6; RCS = 0.1; G = 45;  
%   snr = radareqsnr(lambda,RNG,Pt,tau,'rcs',RCS,'gain',G)
%
%   See also phased, radareqrng, radareqpow.

%   Copyright 2007-2018 The MathWorks, Inc.

%   References
%     [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%         McGraw-Hill, 2001
%     [2] Nicholas J. Willis, Bistatic Radar, Scitech Publishing, 2005


%#codegen
%#ok<*EMCA>
phased.internal.narginchk(4,12,nargin);
[RCS, Ts, Gain, Loss] = parseInput(lambda,RNG,Pt,tau,varargin{:}); % Validate required inputs.

if numel(Gain)==1 % monostatic
    Gtx = db2pow(Gain);
    Grx = Gtx;
else               % bistatic
    % For bistatic radar, split up the Tx and Rx gains.
    Gain = db2pow(Gain);
    Gtx = Gain(1);
    Grx = Gain(2);
end

if numel(RNG)==1  % monostatic
    Rtx = RNG;
    Rrx = RNG;
else               % bistatic
    % For bistatic radar, split up the Tx and Rx ranges.
    Rtx = RNG(1);
    Rrx = RNG(2);
end
k   = physconst('Boltzmann'); % J/K 
L   = db2pow(Loss);
snr = (Pt*tau*Gtx*RCS*Grx*lambda^2)/((4*pi)^3*k*Ts*L*Rtx^2*Rrx^2); % dB
SNRdb = pow2db(snr);

%-------------------------------------------
function [RCS, Ts, Gain, Loss] = parseInput(lambda,RNG,Pt,tau,varargin)

eml_assert_no_varsize(1:nargin,lambda,RNG,Pt,tau,varargin{:});
funName = 'radareqsnr';

validateattributes(lambda,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'LAMBDA',1);

validateattributes(RNG,{'double'},{'positive','nonempty'},...
    funName,'RANGE',2);

cond = numel(RNG) <= 2;
if ~cond
    coder.internal.assert(cond,'phased:radareq:invalidDimension', 'RANGE');
end

validateattributes(Pt,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'Pt',3);

validateattributes(tau,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'Tau',4);

  % Define default values for optional inputs.
  defaultRCS     = 1;      % m^2
  defaultTs      = 290;    % Kelvin
  defaultGain    = 20;     % dB
  defaultLoss    = 0;      % dB

 if isempty(coder.target)
     p = inputParser;
     p.addParameter('RCS',defaultRCS);   
     p.addParameter('Ts',defaultTs); 
     p.addParameter('Gain',defaultGain); 
     p.addParameter('Loss',defaultLoss);  
     p.parse(varargin{:});
     RCS = p.Results.RCS;
     Ts = p.Results.Ts;
     Gain = p.Results.Gain;
     Loss = p.Results.Loss;     
 else

     parms = struct('RCS',uint32(0), ...
           'Ts',uint32(0), ...
           'Gain',uint32(0), ...
           'Loss',uint32(0));

     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     RCS = eml_get_parameter_value(pstruct.RCS,defaultRCS,varargin{:});
     Ts = eml_get_parameter_value(pstruct.Ts,defaultTs,varargin{:});
     Gain = eml_get_parameter_value(pstruct.Gain,defaultGain,varargin{:});
     Loss = eml_get_parameter_value(pstruct.Loss,defaultLoss,varargin{:});     
 end

radareqvalidateoptionalinput(funName,RCS,Ts,Gain,Loss); % Validate optional inputs.

% [EOF]
