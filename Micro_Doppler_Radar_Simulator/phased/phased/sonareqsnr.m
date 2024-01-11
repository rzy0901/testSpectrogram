function SNRdB = sonareqsnr(SL,NL,DI,TL,varargin) 
%sonareqsnr Estimate SNR using sonar equation
%   SNR = sonareqsnr(SL,NL,DI,TL) returns the signal-to-noise- ratio, SNR
%   in dB. SL is the source level in dB//1uPa. NL is the noise Level in
%   dB//1uPa. DI is the sonar directivity index in dB and TL is the
%   transmission loss in dB.
%     
%   sonareqsnr(SL,NL,DI,TL,TS) specifies the target strength, TS in dB. Use
%   this syntax when modeling a signal echo in an active sonar system. If
%   omitted, the sonar system is assumed to be operating in passive mode.
%
%   SNR,SL,NL,DI,TS are real scalars. TL is a positive scalar.
%   
%   % Examples:
%
%   % Example 1: 
%   %   Estimate the SNR required for a sonar with source level of 220 
%   %   dB//1uPa. Undergoing a transmission loss of 85 dB. The directivity
%   %   index is 20 dB. Assuming a target with target strength of 25 dB. 
%   %   The noise level is 73 dB//1uPa.
%   
%   SL = 220; 
%   NL = 73; 
%   DI = 20; 
%   TL = 85; 
%   TS = 25; 
%   snr = sonareqsnr(SL,NL,DI,TL,TS)
%
%   % Example 2: 
%   %   Estimate the SNR required for a sonar to detect a transmission with
%   %   source level of 195 dB//1uPa and transmission loss of 103 dB. The  
%   %   directivity index is 0 dB. Assuming a noise level of 91 dB//1uPa.
%          
%   SL = 195;
%   NL = 91; 
%   DI = 0; 
%   TL = 103; 
%   snr = sonareqsnr(SL,NL,DI,TL)
%
%   See also phased, sonareqtl, sonareqsl.
 
%   Copyright 2017 The MathWorks, Inc

%   References
%   [1]
%   http://www.dosits.org/science/advancedtopics/sonarequation/sonarexample/
%   [2]
%   http://www.dosits.org/science/advancedtopics/sonarequation/psexample/

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(4,5,nargin);
TS = parseInput(SL,NL,DI,TL,varargin{:}); % Validate required inputs.

if (nargin == 5)                    % 5 inputs in case of  Active   
    SNRdB = SL+TS-2*TL-(NL-DI);   
elseif (nargin == 4)                % 4 inputs in case of Passive  
    SNRdB = SL-TL-(NL-DI);
end

%-------------------------------------------
function TS = parseInput(SL,NL,DI,TL,varargin)

funName = 'sonareqsnr';

validateattributes(SL,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'SL',1);

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
