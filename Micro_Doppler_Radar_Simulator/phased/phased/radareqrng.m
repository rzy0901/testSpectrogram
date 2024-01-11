function rng = radareqrng(lambda,SNRdb,Pt,tau,varargin)
%radareqrng Radar range equation
%   RNG = radareqrng(LAMBDA,SNR,Pt,Tau) estimates the maximum theoretical
%   detectable radar range based on the wavelength, LAMBDA (meters),
%   signal-to-noise-ratio, SNR (dB), peak power, Pt (Watts), and the pulse
%   width, Tau (seconds).
%
%   For bistatic radars, the range returned is the geometric mean of the
%   range from the transmitter to the target and from the target to the
%   receiver. The range returned is RNG = sqrt(TxRng*RxRng).
%
%   radareqrng(...,'RCS',RCS) specifies the target's radar cross section
%   RCS in square meters. If omitted RCS = 1 m^2. The function assumes a
%   non-fluctuating target (Swerling case 0).
%
%   radareqrng(...,'Ts',Ts) specifies the system noise temperature Ts
%   in kelvin (K). If omitted Ts = 290 K.
%
%   radareqrng(...,'Gain',G) specifies the antenna gain G in decibels (dB).
%   If omitted G = 20 dB. For a bistatic radar, specify the gain as a two
%   element row vector to represent the transmit antenna and receive
%   antenna gain, G = [TxGain RxGain].
%
%   radareqrng(...,'Loss',L) specifies system losses L in decibels (dB). If
%   omitted L = 0 dB.
%
%   radareqrng(...,'unitstr',UNITS) specifies the units of the returned 
%   range. UNITS can be one of the following strings: 'km' (kilometer),
%   'm' (meter), 'mi' (statue mile), or 'nmi' (nautical mile in the US). If
%   omitted UNITS = 'm'.
%
%   % Example:
%   %   Estimate the range, in miles, which a radar can detect a target
%   %   with an RCS of 0.1 m^2. The radar operates at a frequency of 5.6 
%   %   GHz and transmits at a peak power Pt of 1.5 MW. It uses rectangular
%   %   waveform with a pulse width of 0.2 microseconds. The antenna gain 
%   %   is 45 dB and the minimum detectable SNR is 20 dB. Assume no losses.
%
%   lambda = 3e8/5.6e9; SNR = 20; Pt = 1.5e6; 
%   tau = 0.2e-6; RCS = 0.1; G = 45;  
%   radar_rng = radareqrng(lambda,SNR,Pt,tau,'rcs',RCS,'gain',G,...
%                   'unitstr','mi')
%
%   See also phased, radareqsnr, radareqpow.

%   Copyright 2007-2018 The MathWorks, Inc.

%   References
%     [1] Merrill Skolnik, Introduction to Radar Systems, 3rd Ed.,
%         McGraw-Hill, 2001
%     [2] Nicholas J. Willis, Bistatic Radar, Scitech Publishing, 2005

%#codegen
%#ok<*EMCA>
phased.internal.narginchk(4,14,nargin);
[RCS, Ts, Gain, Loss, unitStr] = parseInput(lambda,SNRdb,Pt,tau,varargin{:}); % Validate required inputs.

if numel(Gain)==1 % monostatic
    Gtx = db2pow(Gain);
    Grx = Gtx;
else               % bistatic
    % For bistatic radar, split up the Tx and Rx gains.
    Gain = db2pow(Gain);
    Gtx = Gain(1);
    Grx = Gain(2);
end

L   = db2pow(Loss);
k   = physconst('Boltzmann'); % J/K 
SNR = db2pow(SNRdb);
rng = ( (Pt*tau*Gtx*RCS*Grx*lambda^2)/((4*pi)^3*k*Ts*L*SNR) )^(1/4); % meters

% Convert to desired units.
switch unitStr
    case 'm'   % meter
        % noop
        
    case 'km'  % kilometers 
        rng = rng/1000; 
        
    case 'nmi' % nautical miles
        rng = rng/1000*0.539956803;
        
    case 'mi'  % miles
        rng = rng/1000*0.621371192; 
end

%-------------------------------------------

function [RCS, Ts, Gain, Loss, unitstr] = parseInput(lambda,SNRdb,Pt,tau,varargin)

eml_assert_no_varsize(1:nargin,lambda,SNRdb,Pt,tau,varargin{:});
funName = 'radareqrng';

validateattributes(lambda,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'LAMBDA',1);

validateattributes(SNRdb,{'double'},{'nonempty','scalar','real'},...
    funName,'SNR',2);

validateattributes(Pt,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'Pt',3);

validateattributes(tau,{'double'},{'positive','nonempty','scalar','finite'},...
    funName,'Tau',4);


% Define default values for optional inputs.
defaultRCS     = 1;      % m^2
defaultTs      = 290;    % Kelvin
defaultGain    = 20;     % dB
defaultLoss    = 0;      % dB
defaultunitStr = 'm';    % meters

 if isempty(coder.target)
     p = inputParser;
     p.addParameter('RCS',defaultRCS);   
     p.addParameter('Ts', defaultTs); 
     p.addParameter('Gain',defaultGain); 
     p.addParameter('Loss',defaultLoss);  
     p.addParameter('unitstr',defaultunitStr);       
     p.parse(varargin{:});
     RCS = p.Results.RCS;
     Ts = p.Results.Ts;
     Gain = p.Results.Gain;
     Loss = p.Results.Loss;     
     unitstr = p.Results.unitstr;     
 else

     parms = struct('RCS',uint32(0), ...
           'Ts',uint32(0), ...
           'Gain',uint32(0), ...
           'Loss',uint32(0), ...
           'unitstr',uint32(0));
     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     RCS = eml_get_parameter_value(pstruct.RCS,defaultRCS,varargin{:});
     Ts = eml_get_parameter_value(pstruct.Ts,defaultTs,varargin{:});
     Gain = eml_get_parameter_value(pstruct.Gain,defaultGain,varargin{:});
     Loss = eml_get_parameter_value(pstruct.Loss,defaultLoss,varargin{:});     
     unitstr = eml_get_parameter_value(pstruct.unitstr,defaultunitStr,varargin{:});          
 end

radareqvalidateoptionalinput(funName,RCS,Ts,Gain,Loss); % Validate optional inputs.
     
 unitstr = ...
    validatestring(unitstr,{'m','km','nmi','mi'},'radareqrng','unitstr');


% [EOF]

