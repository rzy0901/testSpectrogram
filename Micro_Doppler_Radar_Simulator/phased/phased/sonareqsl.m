function SL = sonareqsl(SNR,NL,DI,TL,varargin)
%sonareqsl Estimates source level using sonar equation
%   SL = sonareqsl(SNR,NL,DI,TL) returns the source level, SL in dB//1uPa.
%   SNR is signal to noise ratio in dB. NL is noise level in dB//1uPa. DI
%   is the sonar directivity index in dB and TL is transmission loss in dB.
%
%   sonareqsl(SNR,NL,DI,TL,TS) specifies the target strength, TS in dB. Use
%   this syntax when modeling a signal echo in an active sonar system. If
%   omitted, the sonar system is assumed to be operating in passive mode.
%
%   SL,SNR,NL,DI,TS are real scalars. TL is a positive scalar.
%     
%   % Examples:
%
%   % Example 1: 
%   %   Estimate the source level necessary for detection by sonar with
%   %   SNR of 22 dB and transmission loss of 85 dB.The directivity index
%   %   is 20 dB. Assuming a target of target strength 25 dB. The noise
%   %   level is 73 dB//1uPa.
%   
%   SNR = 22; 
%   NL = 73;
%   DI = 20;
%   TL = 85; 
%   TS = 25; 
%   sl = sonareqsl(SNR,NL,DI,TL,TS)
%
%   % Example 2:
%   %   Estimate the source level of a transmitter required, so that the
%   %   transmission can be detected by a receiver with SNR of 1 dB.The 
%   %   transmission loss is 103 dB.The directivity index is 0 dB. Assume
%   %   noise level of 91 dB//1uPa.
%       
%   SNR = 1;
%   NL = 91;
%   DI = 0; 
%   TL = 103; 
%   sl = sonareqsl(SNR,NL,DI,TL)
%
%   See also phased, sonareqtl, sonareqsnr.
 
%   Copyright 2017 The MathWorks, Inc.
 
%   References
%   [1]
%   http://www.dosits.org/science/advancedtopics/sonarequation/sonarexample/
%   [2]
%   http://www.dosits.org/science/advancedtopics/sonarequation/psexample/

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(4,5,nargin);
TS = parseInput(SNR,NL,DI,TL,varargin{:}); % Validate required inputs.

if (nargin == 5)                            % 5 inputs in case of Active    
    SL = SNR-TS+2*TL+(NL-DI);    
elseif (nargin == 4)                        % 4 inputs in case of Passive    
    SL = SNR+TL+(NL-DI);
end

%-------------------------------------------
function TS = parseInput(SNR,NL,DI,TL,varargin)

funName = 'sonareqsl';

validateattributes(SNR,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'SNR',1);

validateattributes(NL,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'NL',2);

validateattributes(DI,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'DI',3);

validateattributes(TL,{'double'},{'positive','nonempty','scalar','real','finite'},...
    funName,'TL',4);

if ~isempty(varargin)
    TS = varargin{1};
    validateattributes(TS,{'double'},{'nonempty','scalar','real','finite'},funName,'TS',5);
else
    TS = 0;
end

%[EOF]
