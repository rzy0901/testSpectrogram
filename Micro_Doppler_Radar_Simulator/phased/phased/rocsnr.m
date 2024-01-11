function varargout =rocsnr(SNR, varargin)
%rocsnr Receiver operating characteristic curves on varying SNR
%   [Pd, Pfa] = rocsnr(SNRdb) returns the calculated receiver operating
%   characteristic (ROC) curve, i.e., probability of detection (Pd) vs.
%   probability of false alarm (Pfa), for a given single pulse signal to
%   noise ratio (SNRdb) specified in dB. Pd and Pfa are returned in
%   columns. 
%
%   If a vector of SNR values are specified, Pd will be a matrix with each
%   column corresponding to a given SNRdb value. Pfa is always a column
%   vector and is shared by all SNRdb values.
%
%   [Pd, Pfa] = rocsnr(...,'SignalType',TYPE) specifies the type of the
%   received signal when calculating the ROC curve. TYPE can be one of the
%   following: 'Real' | 'NonfluctuatingCoherent' |
%   'NonfluctuatingNoncoherent' | 'Swerling1' | 'Swerling2' | 'Swerling3' |
%   'Swerling4'. The default TYPE is 'NonfluctuatingCoherent'. The noise is
%   assumed to be white Gaussian.
%
%   [Pd, Pfa] = rocsnr(...,'MinPfa',MINPFA) specifies the minimum Pfa
%   included in the ROC calculation. The default value of MINPFA is 1e-10.
%
%   [Pd, Pfa] = rocsnr(...,'MaxPfa',MAXPFA) specifies the maximum Pfa
%   included in the ROC calculation. The default value of MAXPFA is 1.
%
%   [Pd, Pfa] = rocsnr(...,'NumPoints',NUMPTS) specifies the number of
%   points, i.e., steps, used to calculate the ROC curves. Note that the
%   points are evenly spaced in the log scale, i.e., in the exponent of
%   Pfa. The default value of NUMPTS is 101.
%
%   [Pd, Pfa] = rocsnr(...,'NumPulses',NUMINT) integrates NUMINT pulses in
%   calculating the ROC curves. The default value of NUMINT is 1, i.e., no
%   pulse integration.
%
%   rocsnr(...) plots the receiver operating characteristic curves.
%
%   % Example:
%   %   Plot the ROC curve for the nonfluctuating-coherent signal type
%   %   using a set of given SNR values.
%
%   d = [0 3 10 13]; % SNR in dB
%   rocsnr(d,'SignalType','NonfluctuatingCoherent');
%
%   See also phased, rocpfa, npwgnthresh.

%   Copyright 2010 The MathWorks, Inc.

%   Reference 
%   [1] Mark Richards, Fundamentals of Radar Signal Processing
%   [2] Robert McDonough and Anthony Whalen, Detection of Signals in Noise

%#codegen
%#ok<*EMCA>

cond = nargin > 0;
if ~cond
    coder.internal.assert(cond,'MATLAB:narginchk:notEnoughInputs');
end
cond = nargout < 3;
if ~cond
    coder.internal.assert(cond,'MATLAB:nargoutchk:tooManyOutputs');
end

[SignalType, MinPfa, MaxPfa, NumPoints, NumPulses] = ...
    parseInput(SNR,varargin{:});


d = db2pow(SNR);
MinPfa = log10(MinPfa);  % convert to log scale for even space sampling
MaxPfa = log10(MaxPfa);
Pfa = logspace(MinPfa,MaxPfa,NumPoints);

Pd = privrocpdcalc(Pfa,d,NumPulses,SignalType);

if ~nargout
    privrocplot('rocsnr',Pfa,Pd,SNR,@semilogx,'P_{fa}','SNR','dB',SignalType);
else
    varargout{1} = Pd;
end
if nargout > 1
    varargout{2} = Pfa(:);
end

%--------------------------------------------------------------------------

function [SignalType, MinPfa, MaxPfa, NumPoints, NumPulses] = ...
    parseInput(SNR,varargin)

eml_assert_no_varsize(1:nargin, SNR, varargin{:});
% Define default values for optional inputs.
defaultSignalType = 'NonfluctuatingCoherent';
defaultMinPfa = 1e-10;
defaultMaxPfa = 1;
defaultNumPoints = 101;
defaultNumPulses = 1;

 if isempty(coder.target)
     p = inputParser;
     p.addParameter('SignalType',defaultSignalType);   
     p.addParameter('MinPfa',defaultMinPfa); 
     p.addParameter('MaxPfa',defaultMaxPfa); 
     p.addParameter('NumPoints',defaultNumPoints);  
     p.addParameter('NumPulses',defaultNumPulses);  
     p.parse(varargin{:});
     SignalType = p.Results.SignalType;
     MinPfa = p.Results.MinPfa;
     MaxPfa = p.Results.MaxPfa;
     NumPoints = p.Results.NumPoints;     
     NumPulses = p.Results.NumPulses;          
 else

     parms = struct('SignalType',uint32(0), ...
           'MinPfa',uint32(0), ...
           'MaxPfa',uint32(0), ...
           'NumPoints',uint32(0), ...                    
           'NumPulses',uint32(0));

     pstruct = eml_parse_parameter_inputs(parms,[],varargin{:});
     SignalType = ...
         eml_get_parameter_value(pstruct.SignalType,defaultSignalType,varargin{:});
     MinPfa = eml_get_parameter_value(pstruct.MinPfa,defaultMinPfa,varargin{:});
     MaxPfa = eml_get_parameter_value(pstruct.MaxPfa,defaultMaxPfa,varargin{:});
     NumPoints = eml_get_parameter_value(pstruct.NumPoints,defaultNumPoints,varargin{:});     
     NumPulses = eml_get_parameter_value(pstruct.NumPulses,defaultNumPulses,varargin{:});          
 end


SignalType = validatestring(SignalType,{'Real','NonfluctuatingCoherent',...
    'NonfluctuatingNoncoherent','Swerling1','Swerling2','Swerling3',...
    'Swerling4'},'rocsnr','SignalType');


 validateattributes(SNR,{'numeric'},{'vector','real'},...
     'rocsnr','SNR');
 
 validateattributes(MinPfa,{'numeric'},{'scalar','real','>',0,'<',1},...
     'rocsnr','MinPfa');
 
 validateattributes(MaxPfa,{'numeric'},{'scalar','real','>',0,'<=',1},...
     'rocsnr','MaxPfa');
 
 cond = MinPfa < MaxPfa;
 if ~cond
     coder.internal.assert(cond,'phased:rocsnr:InvalidInput');
 end
 
 validateattributes(NumPoints,{'numeric'},{'scalar','positive','integer'},...
     'rocsnr','NumPoints');
 
 validateattributes(NumPulses,{'numeric'},{'scalar','positive','integer'},...
     'rocsnr','NumPulses');



% [EOF]
