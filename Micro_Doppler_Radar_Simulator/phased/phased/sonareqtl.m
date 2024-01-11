function TL = sonareqtl(SL,SNR,NL,DI,varargin)  
%sonareqtl Estimate Transmission loss using sonar equation
%   TL = sonareqtl(SL,SNR,NL,DI) returns the transmission loss, TL in dB.
%   SL is the source level in dB//1uPa. SNR is the signal-to-noise ratio in
%   dB. NL is the noise level in dB//1uPa and DI is the sonar directivity
%   index in dB.
%    
%   sonareqtl(SL,SNR,NL,DI,TS) specifies the target strength, TS in dB. Use
%   this syntax when modeling a signal echo in an active sonar system. If
%   omitted, the sonar system is assumed to be operating in passive mode.
%   
%   SL,SNR,NL,DI,TS are real scalars. TL is a positive scalar.
%
%   % Examples:
%
%   % Example 1: 
%   %   Estimate the one way tranmsission loss for a sonar with source
%   %   level 220 dB//1uPa and SNR of 22 dB. The directivity index is 20 dB
%   %   dB. Assume a target strength of 25 dB and noise level of 73 
%   %   dB//1uPa.
%          
%   SL = 220;
%   SNR = 22; 
%   NL = 73; 
%   DI = 20; 
%   TS = 25;
%   tl = sonareqtl(SL,SNR,NL,DI,TS)
%
%   % Example 2: 
%   %   Estimate the transmission loss for a passive sonar with SNR 1 dB 
%   %   to detect a target with source level of 195 dB//1uPa of 1 dB. The
%   %   directivity index is 0 dB. Assuming a noise level of 91 dB//1uPa.  
%          
%   SL = 195; 
%   SNR = 1; 
%   NL = 91; 
%   DI = 0;
%   tl = sonareqtl(SL,SNR,NL,DI)
%
%   See also phased, sonareqsl, sonareqsnr.
 
%   Copyright 2017 The MathWorks, Inc.
 
%   References
%   [1]
%   http://www.dosits.org/science/advancedtopics/sonarequation/sonarexample/
%   [2]
%   http://www.dosits.org/science/advancedtopics/sonarequation/psexample/

%#codegen
%#ok<*EMCA>

phased.internal.narginchk(4,5,nargin);
TS = parseInput(SL,SNR,NL,DI,varargin{:}); % Validate required inputs.

if (nargin == 5)               % 5 inputs in case of  Active  
    TL = (SL+TS-SNR-(NL-DI))/2;  
elseif (nargin == 4)          % 4 inputs in case of Passive
    TL = SL-SNR-(NL-DI);
end

if ~(TL>0)                  % Validate the output
    error( 'MATLAB:sonareqtl:expectedPositive','Transmission loss should be nonzero and positive,check inputs');
end

%-------------------------------------------
function TS = parseInput(SL,SNR,NL,DI,varargin)

funName = 'sonareqtl';

validateattributes(SL,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'SL',1);

validateattributes(SNR,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'SNR',2);

validateattributes(NL,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'NL',3);

validateattributes(DI,{'double'},{'nonempty','scalar','real','finite'},...
    funName,'DI',4);

if ~isempty(varargin)
    TS = varargin{1};
    validateattributes(TS,{'double'},{'nonempty','scalar','finite'},funName,'TS',5);
else
    TS = 0;
end

%[EOF]
