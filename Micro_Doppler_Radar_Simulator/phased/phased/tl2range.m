function RNG = tl2range(TL,FREQ,DEPTH)
% tl2range Computes range from transmission loss
%   RNG = tl2range(TL,FREQ,DEPTH) returns range, RNG in meters
%   corresponding to the transmission loss, TL in dB. FREQ is the signal
%   frequency in Hz and DEPTH is the water depth in meters.
%
%   TL, FREQ, DEPTH are positive scalars.
%
%   This function calculates the attenuation coefficient for transmission
%   loss using standard values for salinity (35 ppt), pH (8) and
%   temperature (10 C). The model depth is assumed to be DEPTH/2.
%
%   % Example: 
%   %   Given a sonar operating at frequency of 1000 Hz and transmission
%   %   loss of 58.11 dB, estimate the target range, assuming a water depth
%   %   of 2620 m.
%      
%   TL = 58.11; 
%   freq = 1000; 
%   depth = 2620;
%   rng = tl2range(TL,freq,depth)
%
%   See also phased, range2tl.

%   Copyright 2017 The MathWorks, Inc.

%   References 
%   [1] 
%   http://oalib.hlsresearch.com/PocketBook%203rd%20ed.pdf
%   [2]
%   http://calhoun.nps.edu/bitstream/handle/10945/28967/introductiontoso00copp.pdf?sequence=1

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(3,3,nargin);
parseInput(TL,FREQ,DEPTH); % Validate required inputs.

% To calculate attenuation constant
freqkHz = FREQ*1e-3;    % Convert Frequency in Hz to kHz.
T = 10;                 % Temperature (degree celsius)
S = 35;                 % Salinity (ppt)
pH = 8;                 % pH
z = DEPTH/(2*1000);     % Depth (km)
if (z > 7)
    z = 7;              % Set depth to 7 km if it is greater than 7 km
end

f1 = 0.78*sqrt(S/35)*exp(T/26); % Relaxation frequency for boric acid
f2 = 42*exp(T/17);              % Relaxation frequency for magnesium sulphate

a = 0.106*(f1*(freqkHz)^2/(f1^2+(freqkHz)^2))*exp((pH-8)/0.56)...             % Relaxation of boric acid
    +0.52*(1+(T/43))*(S/35)*(f2*(freqkHz)^2/(f2^2+(freqkHz)^2))*exp(-z/6) ...% Relaxation of magnesium sulphate
    + 4.9*1e-04*exp(-((T/27)+(z/17)))*(freqkHz)^2;              % Shear and bulk viscosity of pure water

RNG = 0;                % Initialization for codegen
TRANSDEPTH = DEPTH/2 ;  % Transition depth, after which spherical spreading becomes cylindrical. 

% Find reference TL based on Transition depth
TLtrans = 20*log10(TRANSDEPTH)+a*TRANSDEPTH*1e-3;

if(TL < TLtrans)          % For spherical spreading loss
    RNG = fzero(@(r)20*log10(r)+a*r*(1e-3)-TL,[0.001 TL/(a*1e-3)]);
elseif (TL >= TLtrans)    % For cylindrical spreading loss
    RNG = fzero(@(r)20*log10(TRANSDEPTH)+10*log10(r/TRANSDEPTH)+....
        a*r*(1e-3)-TL,[0.001 (TL-20*log10(TRANSDEPTH))/(a*1e-3)]);
end

%-----------------------------------
function parseInput(TL,FREQ,DEPTH)

funName = 'tl2range';

validateattributes(TL,{'double'},{'positive','nonempty','real','finite'},funName,'TL',1);
validateattributes(FREQ,{'double'},{'positive','nonempty','real','finite','<=',2000000},funName,'Freq',2);
validateattributes(DEPTH,{'double'},{'positive','nonempty','real','finite','>=',2},...
   funName,'depth',3);

%[EOF]
