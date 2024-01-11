function TL = range2tl(RNG,FREQ,DEPTH)
%range2tl Computes transmission loss from range
%   TL = range2tl(RNG,FREQ,DEPTH) returns the transmission loss, TL in dB
%   from range specified in RNG in meters. FREQ is the signal frequency in
%   Hz and DEPTH is the waterdepth in meters.
%
%   RNG, FREQ, DEPTH are positive scalars.
%
%   This functions calculates the attenuation coefficient for transmission
%   loss using standard values for salinity (35 ppt), pH (8) and
%   temperature (10 C). The model depth is assumed to be DEPTH/2.
%
%   % Example: 
%   %   Given a sonar operating at frequency of 1000 Hz. Estimate the  
%   %   transmission loss for a target range of 800 m, assuming the water 
%   %   depth to be 2620 m.
%
%   rng = 800; 
%   freq = 1000; 
%   depth = 2620;
%   tl = range2tl(rng,freq,depth)
%
%   See also phased, tl2range.

%   Copyright 2017 The MathWorks, Inc.

%   References
%   [1] 
%   http://oalib.hlsresearch.com/PocketBook%203rd%20ed.pdf
%   [2] 
%   http://calhoun.nps.edu/bitstream/handle/10945/28967/introductiontoso00copp.pdfsequence=1

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(3,3,nargin);

parseInput(RNG,FREQ,DEPTH); % Validate required inputs.

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
    + 4.9*1e-04*exp(-((T/27)+(z/17)))*(freqkHz)^2;                % Shear and bulk viscosity of pure water

TL = 0;                 % Initialization for codegen
TRANSDEPTH = DEPTH/2 ;  % Transition depth, after which spherical spreading becomes cylindrical. 

if (RNG <= TRANSDEPTH)        % For spherical spreading loss
    TL = pow2db(RNG^2)+a*RNG*1e-3;
elseif (RNG > TRANSDEPTH)     % For cylindrical spreading loss
    TL = pow2db(RNG) + pow2db(TRANSDEPTH) + a*RNG*1e-3;
end

%-----------------------------------
function parseInput(RNG,FREQ,DEPTH)

funName = 'range2tl';

validateattributes(RNG,{'double'},{'positive','nonempty','real','finite','>=',1},funName,'Rng',1);
validateattributes(FREQ,{'double'},{'positive','nonempty','real','finite','<=',2000000},funName,'Freq',2);
validateattributes(DEPTH,{'double'},{'positive','nonempty','real','finite','>=',2},....
   funName,'transdepth',3);

%[EOF]
